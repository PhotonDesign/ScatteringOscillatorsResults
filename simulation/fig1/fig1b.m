%% *** fig1b
% Escat from elliptical scatterer under plane wave incidence

clc
close all

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
nl = 101; % number of frequencies, change to 2001 for accurate result
[s,weight_s] = legpts(nl,[0.02,0.98]);
wmax = 2*pi*1;
ww = wmax.*s; weight = weight_s.'.*wmax;

chi_m = 4;
radius = 0.8;
% *** specify permittivity distribution of scatter using tanh
c1 = 30;
Chi_xy = @(x,y) chi_m *( ( tanh(c1*(-sqrt((x-ay/2).^2+((y-ay/2)./1.5).^2)+radius)) )+1 )./2;
Chi = Chi_xy(X,Y);

% *** points to query Escat
ptx = [0.3,0.1,0.3,0.6,0.9,1.2,1.5]; pty = [0.5,2.4,2.8,2.8,2.8,2.8,2.9];
npt = length(ptx);
ptindex = nx*(floor(pty/hx)+1) + round(ptx/hx)+1;

% Trapezoidal matrix
B = trapez_mat(ax, ay, nx, ny);

% *** find the grids where the scatter is defined
xscatter = zeros(ny,nx);
yscatter = zeros(ny,nx);
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
[~,~,chi_s] = find(chi_scatter);


es_iw = zeros(nl,npt);
for iw = 1:nl
    tic
    wsim = ww(iw); wl = 2*pi*1/wsim;
    k = 2*pi/wl;
    
    % incident field
    phi = exp(-1j*k*ax/2);
    Einc_xy = @(x,y) exp(1j*k*x) * phi;
    Einc = Einc_xy(X,Y);
    
    % ***
    % parameters for BICGSTAB
    tol = 1e-10;
    maxit = 8000;
    
    quad_order = 6;
    % *** D0 and D1
    [D0_t,D1_t]= planck_win_get_D(k,hx);
    D0=D0_t*hx^2; D1=D1_t*hx^4;
    
    % *** create the GF convolutor
    x2 = -ax:hx:ax;
    y2 = -ay:hy:ay;
    [X2, Y2] = meshgrid(x2,y2);
    r = sqrt(X2.^2 + Y2.^2);
    g = 1j*k^2 / 4 * besselh(0,k*r);
    g(ny,nx) = 0; % set singularity (at the origin) to zero
    
    % *** high-order quadrature correction
    Tao = zeros(size(g));
    switch quad_order
        case 0 % no correction, do nothing
        case 4 % forth order correction
            Tao(ny,nx) = D0;
        case 6 % sixth order correction
            Tao(ny,nx) = D0 - 2*D1/hx^2;
            Tao(ny-1,nx) = 1/2 * D1 / hx^2;
            Tao(ny+1,nx) = 1/2 * D1 / hx^2;
            Tao(ny,nx-1) = 1/2 * D1 / hx^2;
            Tao(ny,nx+1) = 1/2 * D1 / hx^2;
        case 10
            Tao(ny,nx) = D0 - 49/18*D1/hx^2 + 7/9*D2/hx^4 - 1/18*D4/hx^6 - 1/2*D5/hx^6;
            [Tao(ny-1,nx), Tao(ny+1,nx), Tao(ny,nx-1), Tao(ny,nx+1)] = deal( 3/4*D1/hx^2 - 13/48*D2/hx^4 - 19/24*D3/hx^4 + 1/48*D4/hx^6 + 7/24*D5/hx^6 );
            [Tao(ny-2,nx), Tao(ny+2,nx), Tao(ny,nx-2), Tao(ny,nx+2)] = deal( -3/40*D1/hx^2 + 1/12*D2/hx^4 + 1/24*D3/hx^4 -1/120*D4/hx^6 - 1/24*D5/hx^6 );
            [Tao(ny-3,nx), Tao(ny+3,nx), Tao(ny,nx-3), Tao(ny,nx+3)] = deal( 1/180*D1/hx^2 - 1/144*D2/hx^4 + 1/720*D4/hx^6 );
            [Tao(ny-1,nx-1), Tao(ny-1,nx+1), Tao(ny+1,nx-1), Tao(ny+1,nx+1)] = deal( 5/12*D3/hx^4 - 1/6*D5/hx^6 );
    end
    g = g + 1j*k^2/4/hx^2 * Tao;
    
    % iterative method
    afun = @(escat) escat - BTTB_matvec(g, B .* Chi .* reshape(escat, ny, nx) * hx * hy);
    yval = BTTB_matvec(g, B .* Chi .* Einc * hx * hy);
    escat = bicgstab(afun, yval, tol, maxit);
    
    % other physical quantities
    Escat = reshape(escat, ny, nx);
    E = Escat + Einc;
    P = Chi .* E;
    
    %figure; imagesc(x,y,transpose(real(Escat))); colorbar; axis equal; axis tight; set(gca,'YDir','Normal');
    
    es_iw(iw,:) = escat(ptindex).';
    toc
end

save('Escat_ellipse.mat','ww','es_iw')

c6 = [25,159,118]./255;
cblue = 0.9.*[65,105,225]./255;
cyellow = 0.9.*[0.9766, 0.7344, 0.2148];
cpurple = 1.3.*[86, 34, 112]./255;
cgreen = [0, 0.2578, 0];
cred = 1.*[0.7695, 0, 0];


lw = 3.2; fs = 60;
figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.15 0.19 0.7 0.7]);
hold on
%xlim([0,1])
ylim([-3,3])
ylabel('Im T','FontSize',fs+3,'interpreter','tex');
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on');
set(gcf,'position',[40,200,1150,1000]); 
hold on
e1re=plot(ww./2./pi,real(es_iw(:,1)),'linewidth',lw+0.8,'color',cblue);
e1im=plot(ww./2./pi,imag(es_iw(:,1)),'--','linewidth',lw+0.8,'color',cblue);
e3re=plot(ww./2./pi,real(es_iw(:,5)),'linewidth',lw+0.8,'color',cred);
e3im=plot(ww./2./pi,imag(es_iw(:,5)),'--','linewidth',lw+0.8,'color',cred);
e6re=plot(ww./2./pi,real(es_iw(:,2)),'linewidth',lw+0.8,'color',cyellow);
e6im=plot(ww./2./pi,imag(es_iw(:,2)),'--','linewidth',lw+0.8,'color',cyellow);

ylabel('Scattered field, E_{scat}','FontSize',fs+3,'interpreter','tex');
xlabel('Normalized frequency, \omega (c/a)','fontsize',fs,'interpreter','tex','fontname','Helvetica');

box(axes2,'off'); axis on
xticks([0,0.2,0.4,0.6,0.8,1])

set(gcf,'renderer','painters')
print -depsc2 Escat.eps
