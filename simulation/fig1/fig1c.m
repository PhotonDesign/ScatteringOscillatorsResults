%% *** fig1c
% several T matrix elements across frequencies

clc
close all

% *** plotting specs
lw = 2.1; fs = 56;
mmap = zeros(64,3); mmap(:,1) = linspace(1,0.5,64); mmap(:,2) = mmap(:,1); mmap(:,3) = mmap(:,1);

% *** region geo
ax = 3;
ay = ax;
nx_list = [32,64,128,256,512]; ny_list = nx_list;

i = 2;
nx = nx_list(i);
ny = ny_list(i);

% *** grid
hx = ax / (nx - 1); % hx = kh_list(i)/k;
hy = ay / (ny - 1);
x = (0:hx:ax).';
y = (0:hy:ay).';
[X, Y] = meshgrid(x,y);

% *** legpts specs, wavelengths and material
nl = 11; % number of frequencies, change to 2001 for accurate result
[s,weight_s] = legpts(nl,[0.02,0.98]);
wmax = 2*pi*1;
ww = wmax.*s; weight = weight_s.'.*wmax;

chi_m = 4;
radius = 0.8;
% *** specify permittivity distribution of scatter using tanh
c1 = 30;
Chi_xy = @(x,y) chi_m *( ( tanh(c1*(-sqrt((x-ay/2).^2+((y-ay/2)./1.5).^2)+radius)) )+1 )./2;


Chi = Chi_xy(X,Y);
figure; imagesc(x,y,transpose(real(Chi))); colorbar; axis equal; axis tight; set(gca,'YDir','Normal');
colormap(mmap); c = colorbar; c.LineWidth = lw-1; c.Limits = [0 chi_m]; c.FontSize = fs;

start = tic; 
% Trapezoidal matrix
B = trapez_mat(ax, ay, nx, ny);

% *** find the grids where the scatter is defined
xscatter = zeros(ny,nx);
yscatter = zeros(ny,nx);
Einc_scatter = zeros(ny,nx);
chi_scatter = zeros(ny,nx);
chi_filter = 0.02;
ng = nnz(Chi>chi_filter);
scatter_index = zeros(2,ng);
iis = 0;
for iix = 1:nx
    for iiy = 1:ny
        if Chi(iiy,iix) > chi_filter
            xscatter(iiy,iix) = X(iiy,iix);
            yscatter(iiy,iix) = Y(iiy,iix);
            chi_scatter(iiy,iix) = Chi(iiy,iix);
            iis = iis+1;
            scatter_index(:,iis) = [iiy,iix];
        end
    end
end
[~,~,xs] = find(xscatter);
[~,~,ys] = find(yscatter);
[row,col,chi_s] = find(chi_scatter);

%
Teigs = zeros(nl,6);
Tdiags = zeros(nl,ng);

Gdiags = zeros(nl,ng);
invGdiags = zeros(nl,ng);

Toff_12 = zeros(nl,1);
Toff_13 = zeros(nl,1);
Toff_23 = zeros(nl,1);
for iw = 1:nl
    % *** setup
    wsim = ww(iw); wl = 2*pi*1/wsim;

    k = 2*pi/wl;
    tic
    
    % *** build GF on scatter grid
    R = zeros(ng,ng);
    G = zeros(ng,ng);
    for ir = 1:ng
        xx = xs - xs(ir);
        yy = ys - ys(ir);
        rr = sqrt(xx.^2+yy.^2);
        R(ir,:) = rr;
        G(ir,:) = 1j.*k.^2 ./ 4 .* besselh(0,k*rr);
        G(ir,ir) = 0;
    end
    
    quad_order = 6;
    % *** D0 and D1
    [D0_t,D1_t]= planck_win_get_D(k,hx);
    D0=D0_t*hx^2; D1=D1_t*hx^4;
    % *** Duan-Rokhlin correction
    Tao = zeros(ng,ng);
    switch quad_order
        case 0 % no correction, do nothing
        case 4 % forth order correction
            Tao(1:ng+1:end) =  D0.*ones(1,ng);
        case 6 % sixth order correction
            for ir1 = 1:ng
                for ir2 = 1:ng
                    xx = xs(ir1) - xs(ir2);
                    yy = ys(ir1) - ys(ir2);
                    rr = sqrt(xx^2+yy^2);
                    if abs(rr)<1e-14
                        Tao(ir1,ir2) = D0 - 2*D1/hx^2;
                    elseif abs(rr-hx)<1e-14
                        Tao(ir1,ir2) = 1/2 * D1 / hx^2;
                    end
                end
            end
    end
    
    G = G + 1j*k^2/4/hx^2 * Tao;
    
    % *** material matrix
    Xi_vec = -1./chi_s;
    Xi = zeros(ng,ng);
    Xi(1:ng+1:end) = Xi_vec;
    toc
    % *** T matrix
    T = -inv(G*hx^2 + Xi);

    Tdiags(iw,:) = diag(T);
    Teigs(iw,:) = eigs(T,6);
    
    Toff_12(iw) = T(400,700);
    Toff_13(iw) = T(400,1200);
    Toff_23(iw) = T(700,1200);

    
    toc
end

elapsed = toc(start);
save('full_T_offdiag.mat','ww','Tdiags','Toff_12','Toff_13','Toff_23','ng','hx','hy','ax','ay')

%%
testrr = [400,700,800,900,1200];
nr = length(testrr);
xi = zeros(1,nr); yi = zeros(1,nr);
for ir = 1:nr
    testr = testrr(ir);
    Einc_s = zeros(ng,1);
    Einc_s(testr) = 1;
    % convert back to whole-region Einc
    Einc = zeros(ny,nx);
    ind = sub2ind(size(Einc),scatter_index(1,:),scatter_index(2,:));
    Einc(ind) = Einc_s;
    [ri,ci] = find(Einc);
    xi(ir) = (ri-1)*hx - ax/2; 
    yi(ir) = (ci-1)*hy - ay/2;
end

c6 = [25,159,118]./255;
cblue = 0.9.*[65,105,225]./255;
cyellow = 0.9.*[0.9766, 0.7344, 0.2148];
cpurple = 1.3.*[86, 34, 112]./255;
cgreen = [0, 0.2578, 0];
cred = 1.*[0.7695, 0, 0];


%% point1, 2-y-axis, plot just 1 diagonal element, both real and imaginary part

T1 = Tdiags(:,400);

lw = 2.1; fs = 56;
figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.2,1])
yticks([0,0.5])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',cblue);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([3.8,5])
yticks([4,4.5])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',cblue);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');

box(axes2,'on'); axis on

set(gcf,'renderer','painters')
print -depsc2 x1.eps
%% point2, 2-y-axis, plot just 1 diagonal element, both real and imaginary part
close all
T1 = Tdiags(:,700);

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.2,0.5])
yticks([0,0.4])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',cyellow);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([4.0,4.6])
yticks([4.2,4.5])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',cyellow);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 x2.eps

%% point3, 2-y-axis, plot just 1 diagonal element, both real and imaginary part
close all
T1 = Tdiags(:,1200);

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.2,1.4])
yticks([0.5,1])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',cred);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([3.6,5])
yticks([4.0,4.5])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',cred);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 x3.eps


%% off_diag_12, 2-y-axis, both real and imaginary part
close all
T1 = Toff_12;

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.16,0.2])
% yticks([0,0.5])
yticks([-0.1 0 0.1])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',cgreen);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([-0.2,0.16])
%yticks([4,4.5])
yticks([-0.1, 0, 0.1])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',cgreen);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');
%xlabel('Normalized frequency, 1/\lambda','fontsize',fs,'interpreter','tex');

box(axes2,'off'); axis on

set(gcf,'renderer','painters')
print -depsc2 off_12.eps
%% off_diag_13, 2-y-axis, both real and imaginary part
close all
T1 = Toff_13;

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.16,0.16])
yticks([-0.1, 0, 0.1])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',cpurple);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([-0.16,0.16])
% yticks([4,4.5])
yticks([-0.1, 0, 0.1])
xticks([0,0.2,0.4,0.6,0.8,1])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',cpurple);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');
%xlabel('Normalized frequency, 1/\lambda','fontsize',fs,'interpreter','tex');
%xlabel('Normalized frequency, \omega (c/a)','fontsize',fs+4,'interpreter','tex','fontname','Helvetica');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 off_13.eps

%% off_diag_23, 2-y-axis, both real and imaginary part
close all
T1 = Toff_23;
c23 = [0,180,180]./255;

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.2,0.1])
yticks([-0.1, 0])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
d1=plot(ww./2./pi,imag(T1),'--','linewidth',lw+0.8,'color',c23);
set(gca,'fontsize',fs)
hold on

axes1 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
axes1.YAxisLocation = 'right';
axes1.Color = 'none';
set(axes1,'fontsize',fs,'XAxisLocation','top')
hold(axes1,'on');
set(gcf,'position',[40,200,1090,335]); 
hold on
ylim([-0.05, 0.25])
yticks([0, 0.2])
xticks([0,0.2,0.4,0.6,0.8,1])
d1re=plot(ww./2./pi,real(T1),'linewidth',lw+0.8,'color',c23);

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 off_23.eps

