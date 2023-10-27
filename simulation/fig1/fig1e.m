%% *** fig1e
% verify sum rule and show monotonicity theorem for T matrix
% 1. to verify low-frequency sum rule, we perform matrix integration of ImT/omega
% 2. to show monotonicity theorem, we test 2 other geometries bounded by the
% elliptical domain, compute their sum-rule results and compare with that
% of the elliptical scatterer. 


clc
close all

format long

% *** region geo
geo.ax = 3;
geo.ay = geo.ax;
nx_list = [16,64,128,256,512,1024]; ny_list = nx_list;

i = 2; % i = 6 for plotting high-res picture of the scatterer shape
geo.nx = nx_list(i);
geo.ny = ny_list(i);

% *** grid
geo.hx = geo.ax / (geo.nx - 1); 
geo.hy = geo.ay / (geo.ny - 1);
x = (0:geo.hx:geo.ax).';
y = (0:geo.hy:geo.ay).';
[geo.X, geo.Y] = meshgrid(x,y);

% *** legpts specs, wavelengths and material
nl = 101; % number of frequencies, to verify sum rule, need 2001 for accurate result
[s,weight_s] = legpts(nl,[0,2]);
wmax = 2*pi*4; 
ww = wmax.*s; weight = weight_s.'.*wmax;

geo.chi_m = 4;
wp = 2*pi*2; g = wp*0.001;
w0 = wp./2; 
chifunc = @(ww) wp.^2./(w0.^2 - ww.^2 - 1i.*g.*ww);
radius = 0.8;
% *** specify permittivity distribution of scatter using tanh
c1 = 30;

% below are 3 different scatterer geometry functions used in the paper

% 1st geometry in fig1e
Chi_xy = @(x,y) geo.chi_m *( ( tanh(c1*(-sqrt((x-geo.ay/2).^2+((y-geo.ay/2)./1.5).^2)+radius)) )+1 )./2;

% 2nd geometry in fig1e
%Chi_xy = @(x,y) geo.chi_m.*( ((x-geo.ax/2) < (-2/3.*(y-geo.ay/2) + 0.8)).* ((x-geo.ax/2) < (1.2.*(y-geo.ay/2) + 0.8)).* ((x-geo.ax/2) > (2/11.*(y-geo.ay/2) - 2.4/11)) );

% 3rd geometry in fig1e
%Chi_xy = @(x,y) geo.chi_m.*( ((y-geo.ay/2) > (0.2.*sin(2./0.8.*pi.*(x-geo.ax/2)).*sin(0.9.*pi.*(x-geo.ax/2)))) .* ((y-geo.ay/2) < sqrt( 1.44.*(1 - (x-geo.ax/2).^2./0.64)))  ) .* ((x-geo.ax/2)<0.8).* ((x-geo.ax/2)>(-0.8));

geo.Chi = Chi_xy(geo.X,geo.Y);

% *** find the grids where the scatter is defined
chi_filter = 0.8;
geo.ng = nnz(abs(real(geo.Chi))>chi_filter*abs(real(geo.chi_m)));
[geo.xs,geo.ys,geo.chi_s] = find_scatterer(geo,chi_filter);
chi_s_0 = geo.chi_s;

% *** find T matrix for this frequency
wsim = ww(1);
T = compute_T(wsim, geo);


% *** get the coordinates for the ellipse
geo_ellipse.nx = geo.nx; geo_ellipse.ny = geo.ny;
geo_ellipse.ax = geo.ax; geo_ellipse.ay = geo.ay;
geo_ellipse.X = geo.X; geo_ellipse.Y = geo.Y;
geo_ellipse.chi_m = 4;
Chi_ellipse_xy = @(x,y) geo_ellipse.chi_m *( ( tanh(c1*(-sqrt((x-geo.ay/2).^2+((y-geo.ay/2)./1.5).^2)+radius)) )+1 )./2;
geo_ellipse.Chi = Chi_ellipse_xy(geo.X,geo.Y);
chi_filter = 0.8;
geo_ellipse.ng = nnz(abs(real(geo_ellipse.Chi))>chi_filter*abs(real(geo_ellipse.chi_m)));
[geo_ellipse.xs,geo_ellipse.ys,geo_ellipse.chi_s] = find_scatterer(geo_ellipse,chi_filter);


% compute the low-frequency sum-rule integral

sum_rule_integral = zeros(geo.ng,geo.ng);

for iw = 1:nl

    wsim = ww(iw); wl = 2*pi*1/wsim;
    chi_w = chifunc(wsim);
    geo.ch_s = chi_s_0./geo.chi_m.*chi_w;
    T = compute_T(wsim, geo);
    sum_rule_integral = sum_rule_integral + imag(T)./wsim.*weight(iw);

end

% embed the freeform structure T matrix into the ellipse T matrix

T_small = sum_rule_integral;
coord_ellipse = geo_ellipse.xs.*1e5 + geo_ellipse.ys;
coord = geo.xs.*1e5 + geo.ys;

embed_index = zeros(geo.ng,1);

for ii = 1:geo.ng
    if find(coord_ellipse == coord(ii))
        embed_index(ii) = find(coord_ellipse == coord(ii));
    else
        embed_index(ii) = 1;
    end
end

T_big = zeros(geo_ellipse.ng, geo_ellipse.ng);
for ii1 = 1:geo.ng
    for ii2 = 1:geo.ng
        T_big(embed_index(ii1),embed_index(ii2)) = T_small(ii1,ii2);
    end
end

save('freq_wmax4.mat')

%% plot the matrix

mred = [0.75,0,0]; mblue = [0,0,0.75]; mwhite = [1,1,1];
mmap = zeros(64,3);
mmap(1:32,1:3) = linspace(1,0,32).'*mblue + linspace(0,1,32).'*mwhite;
mmap(33:64,1:3) = linspace(1,0,32).'*mwhite + linspace(0,1,32).'*mred;

lw = 2.1; fs = 44; ms = 14;
figure1 = figure;

axes2 = axes('Parent',figure1,'Position',[0.15 0.19 0.7 0.7]);
set(gcf,'color','w');
hold on
xlim([0,3])
ylim([0,3])

T_plot = T_big;
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(T_plot)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;

lim = max(max(abs(T_big)));
caxis([-lim lim])

colormap(mmap);
c = colorbar; c.LineWidth = lw-2; %c.Limits = [-0.4 0.4]; c.FontSize = fs;
%ylabel('Im T','FontSize',fs+3,'interpreter','tex');

yticks(linspace(20,80,0));
xticks(linspace(20,80,0));
set(axes2,'YColor','none','XColor','none','fontsize',fs,'box','off');
set(gca,'fontsize',fs)
set(gcf,'position',[10,10,820,800]);
hold on

%% compute error 
T_small = sum_rule_integral;
msize = size(T_small,1);
normfroT = norm(T_small, 'fro');
normfrodiag = norm(diag(T_small), 'fro'); % not the identity matrix, because we used boundary smoothing
error = abs(normfroT - normfrodiag)./normfrodiag;


%%
function [xs,ys,chi_s] = find_scatterer(geo,chi_filter)

xscatter = 1e6.*ones(geo.ny,geo.nx);
yscatter = 1e6.*ones(geo.ny,geo.nx);
chi_scatter = zeros(geo.ny,geo.nx);

scatter_index = zeros(2,geo.ng);
iis = 0;
for iix = 1:geo.nx
    for iiy = 1:geo.ny
        if abs(real(geo.Chi(iiy,iix)))>chi_filter*abs(real(geo.chi_m))
            xscatter(iiy,iix) = geo.X(iiy,iix);
            yscatter(iiy,iix) = geo.Y(iiy,iix);
            chi_scatter(iiy,iix) = geo.Chi(iiy,iix);

            iis = iis + 1;
            scatter_index(:,iis) = [iiy,iix];
        end
    end
end
locx = xscatter<=geo.ax;
xs = xscatter(locx);
locy = yscatter<=geo.ay;
ys = yscatter(locy);
[~,~,chi_s] = find(chi_scatter);


% *** smooth the geometry and plot on finer grid

nxq = 1024; nyq = nxq;
hxq = geo.ax / (nxq - 1); % hx = kh_list(i)/k;
hyq = geo.ay / (nyq - 1);
xq = (0:hxq:geo.ax).';
yq = (0:hyq:geo.ay).';

[Xq,Yq] = meshgrid(xq,yq);
Chi_this_fine = interp2(geo.X,geo.Y,geo.Chi,Xq,Yq);

Chi_this_fine_plot = zeros(nyq,nxq);
for iix = 1:nxq
    for iiy = 1:nyq
        if Chi_this_fine(iiy,iix) > chi_filter*abs(real(geo.chi_m))
        Chi_this_fine_plot(iiy,iix) = Chi_this_fine(iiy,iix);
        end 
    end
end

lw = 2.1; fs = 47; 
mmap = zeros(64,3); mmap(:,1) = linspace(1,0.65,64); mmap(:,2) = mmap(:,1); mmap(:,3) = mmap(:,1);

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.15 0.19 0.7 0.7]);
hold on
xlim([0,3])
ylim([0,3])
yticks(linspace(0,3,0));
xticks(linspace(0,3,0));
imagesc(xq,yq,transpose(real(Chi_this_fine_plot))); axis equal; axis tight; set(gca,'YDir','Normal');
colormap(mmap); %c = colorbar; c.LineWidth = lw-1; c.Limits = [0 chi_m]; c.FontSize = fs;
%ylabel('Im T','FontSize',fs+3,'interpreter','tex');
%set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on');
set(axes2,'YColor','none','XColor','none','fontsize',fs,'box','off');
set(gca,'fontsize',fs)
hold on


xlim([0,3])
ylim([0,3])

box(axes2,'off'); axis on
set(gcf,'position',[40,200,1000,1000]);

end