


clc
close all

sumT = readmatrix('fig1e.txt');
ng = 72;
Tgeo1_w1 = sumT(:,1:ng);
Tgeo1_w2 = sumT(:,ng+1:ng*2);
Tgeo1_inf = sumT(:,2*ng+1:ng*3);

Tgeo2_w1 = sumT(:,3*ng+1:ng*4);
Tgeo2_w2 = sumT(:,4*ng+1:ng*5);
Tgeo2_inf = sumT(:,5*ng+1:ng*6);

Tgeo3_w1 = sumT(:,6*ng+1:ng*7);
Tgeo3_w2 = sumT(:,7*ng+1:ng*8);
Tgeo3_inf = sumT(:,8*ng+1:ng*9);

mred = [0.75,0,0]; mblue = [0,0,0.75]; mwhite = [1,1,1];
mmap = zeros(64,3);
mmap(1:32,1:3) = linspace(1,0,32).'*mblue + linspace(0,1,32).'*mwhite;
mmap(33:64,1:3) = linspace(1,0,32).'*mwhite + linspace(0,1,32).'*mred;

lw = 2.1; fs = 44; ms = 14;


%% plot all panels (not pretty)
figureall = figure;

axes1 = axes('Parent',figureall,'Position',[0.15 0.19 0.7 0.7]);
set(gcf,'color','w');

%
subplot(3,3,1)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo1_w1)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo1_w1)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,4)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo1_w2)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo1_w2)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,7)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo1_inf)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo1_inf)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

%
subplot(3,3,2)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo2_w1)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo2_w1)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,5)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo2_w2)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo2_w2)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,8)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo2_inf)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo2_inf)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

%
subplot(3,3,3)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo3_w1)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo3_w1)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,6)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo3_w2)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo3_w2)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

subplot(3,3,9)
imagesc(1:1:geo_ellipse.ng,1:1:geo_ellipse.ng,flipud(Tgeo3_inf)); axis equal; axis tight; set(gca,'YDir','Normal'); colorbar;
lim = max(max(abs(Tgeo3_inf)));
caxis([-lim lim])
colormap(mmap);
c = colorbar; c.LineWidth = lw-2;

%% individual plot
T_plot = Tgeo1_w1; % choose which matrix to plot, Tgeo1_w1, Tgeo1_w2, Tgeo1_inf, Tgeo2_w1, Tgeo2_w2, Tgeo2_inf, Tgeo3_w1, Tgeo3_w2, Tgeo3_inf

figure1 = figure;

axes2 = axes('Parent',figure1,'Position',[0.15 0.19 0.7 0.7]);
set(gcf,'color','w');
hold on
xlim([0,3])
ylim([0,3])

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