clc
close all

w_T = readmatrix('fig1c.txt');
ww = real(w_T(:,1));
T1 = (w_T(:,2));
T2 = (w_T(:,3));
T3 = (w_T(:,4));
Toff_12 = (w_T(:,5));
Toff_13 = (w_T(:,6));
Toff_23 = (w_T(:,7));

c6 = [25,159,118]./255;
cblue = 0.9.*[65,105,225]./255;
cyellow = 0.9.*[0.9766, 0.7344, 0.2148];
cpurple = 1.3.*[86, 34, 112]./255;
cgreen = [0, 0.2578, 0];
cred = 1.*[0.7695, 0, 0];


%% point1, 2-y-axis, plot just 1 diagonal element, both real and imaginary part


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
d1=plot(ww./2./pi,imag(T2),'--','linewidth',lw+0.8,'color',cyellow);
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
d1re=plot(ww./2./pi,real(T2),'linewidth',lw+0.8,'color',cyellow);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 x2.eps

%% point3, 2-y-axis, plot just 1 diagonal element, both real and imaginary part
close all


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
d1=plot(ww./2./pi,imag(T3),'--','linewidth',lw+0.8,'color',cred);
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
d1re=plot(ww./2./pi,real(T3),'linewidth',lw+0.8,'color',cred);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 x3.eps


%% off_diag_12, 2-y-axis, both real and imaginary part
close all

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
d1=plot(ww./2./pi,imag(Toff_12),'--','linewidth',lw+0.8,'color',cgreen);
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
d1re=plot(ww./2./pi,real(Toff_12),'linewidth',lw+0.8,'color',cgreen);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');
%xlabel('Normalized frequency, 1/\lambda','fontsize',fs,'interpreter','tex');

box(axes2,'off'); axis on

set(gcf,'renderer','painters')
print -depsc2 off_12.eps
%% off_diag_13, 2-y-axis, both real and imaginary part
close all

figure1 = figure;
set(gcf,'color','w');
axes2 = axes('Parent',figure1,'Position',[0.01 0.01 0.99 0.99]);
hold on
xlim([0,1])
ylim([-0.16,0.16])
yticks([-0.1, 0, 0.1])
xticks([0, 0.2, 0.4, 0.6 0.8 1.0])
set(axes2,'YColor','k','XColor','k','fontsize',fs,'box','on','XAxisLocation','bottom');
d1=plot(ww./2./pi,imag(Toff_13),'--','linewidth',lw+0.8,'color',cpurple);
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
d1re=plot(ww./2./pi,real(Toff_13),'linewidth',lw+0.8,'color',cpurple);
%ylabel('Re T','FontSize',fs+3,'interpreter','tex');
%xlabel('Normalized frequency, 1/\lambda','fontsize',fs,'interpreter','tex');
%xlabel('Normalized frequency, \omega (c/a)','fontsize',fs+4,'interpreter','tex','fontname','Helvetica');

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 off_13.eps

%% off_diag_23, 2-y-axis, both real and imaginary part
close all
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
d1=plot(ww./2./pi,imag(Toff_23),'--','linewidth',lw+0.8,'color',c23);
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
d1re=plot(ww./2./pi,real(Toff_23),'linewidth',lw+0.8,'color',c23);

box(axes2,'off'); axis on


set(gcf,'renderer','painters')
print -depsc2 off_23.eps

