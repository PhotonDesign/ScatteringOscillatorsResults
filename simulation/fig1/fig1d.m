%% *** fig1d
% eigenvalues of T matrix across frequencies

clc
close all

% *** region geo
ax = 3;
ay = ax;
nx_list = [16,64,128,256,512,1024]; ny_list = nx_list;

i = 1;
nx = nx_list(i);
ny = ny_list(i);

% *** grid
hx = ax / (nx - 1); % hx = kh_list(i)/k;
hy = ay / (ny - 1);
x = (0:hx:ax).';
y = (0:hy:ay).';
[X, Y] = meshgrid(x,y);

% *** legpts specs, wavelengths and material
nl = 101; % number of frequencies, change to 2001 for accurate result
%[s,weight_s] = legpts(nl,[0.02,0.98]);
[s,weight_s] = legpts(nl,[0.002,0.998]);
wmax = 2*pi*1;
ww = wmax.*s; weight = weight_s.'.*wmax;

chi_m = 4;
radius = 0.8;
% *** specify permittivity distribution of scatter using tanh
c1 = 30;
Chi_xy = @(x,y) chi_m *( ( tanh(c1*(-sqrt((x-ay/2).^2+((y-ay/2)./1.5).^2)+radius)) )+1 )./2;
Chi = Chi_xy(X,Y);

start = tic; 
% Trapezoidal matrix
B = trapez_mat(ax, ay, nx, ny);

% *** find the grids where the scatter is defined
xscatter = zeros(ny,nx);
yscatter = zeros(ny,nx);
Einc_scatter = zeros(ny,nx);
chi_scatter = zeros(ny,nx);
chi_filter = 0.02;
ng = nnz(abs(real(Chi))>chi_filter*abs(real(chi_m)));
scatter_index = zeros(2,ng);
iis = 0;
for iix = 1:nx
    for iiy = 1:ny
        if abs(real(Chi(iiy,iix)))>chi_filter*abs(real(chi_m))
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
ImTeig = zeros(nl,ng);

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

    ImTeig(iw,:) = eig(imag(T));
    
    toc
end

elapsed = toc(start);
save('eig_ImT_ellipse1.mat','ImTeig','ww','ng')



%% plot tall plot
close all
load('eig_ImT_ellipse1.mat')
close all

ImTeig_plot = real(ImTeig);
c1 = 1.3*[0.7695, 0, 0]; 
lw = 2.8; fs = 75;
figure1 = figure;
set(figure1,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
set(gcf,'color','w');


axes2 = axes('Parent',figure1,'Position',[0.19 0.16 0.7 0.8]);
set(axes2,'fontsize',fs)

xticks(linspace(0,1,6));
set(axes2,'YColor','none','XColor','k','fontsize',fs,'box','off');
set(gca, 'TickDir', 'out')


axes1 = axes('Parent',figure1,'Position',[0.19 0.16 0.7 0.80]);
set(axes1,'fontsize',fs)
hold(axes1,'on');
box(axes1,'on');
hold on
for ip = 1:ng
    plot(ww./2./pi,ImTeig_plot(:,ip)./ww.^2.*(2.*pi).^2,'linewidth',lw+1.0,'color',c1.*(log(ng+1-ip)./log(ng+1e5)).*(ip>10) + c1.*log(ng+1-ip)./log(ng+1e1).*(ip<=3)+c1.*log(ng+1-ip)./log(ng+1e3).*(ip>3).*(ip<=10));

end
set(gca, 'YScale', 'linear')
xlim([0.002,0.998])
xlim([0,1])
ylim([-1200,1200])
yticks(linspace(-1200,1200,7));

ylim([-400,800])
yticks(linspace(-1200,1200,7));

xticks(linspace(0,1,0));
set(gca,'xaxislocation','origin','fontsize',fs-5)
ax = gca; ax.XAxis.Exponent = 0;
hold on

ybars = [-400 -1]; cshade = 1.*[0.6953, 0.7930, 0.9648];
patch([min(xlim) max(xlim) max(xlim) min(xlim)],[ybars(1) ybars(1),ybars(2) ybars(2)],cshade,'LineStyle','none');


set(gcf,'position',[10,10,1550,800]);

set(axes1,'YColor','k','XColor','k','fontsize',fs,'box','on');

set(gca, 'TickDir', 'out')
set(gcf,'renderer','painters')
print -depsc2 eigenvalues.eps