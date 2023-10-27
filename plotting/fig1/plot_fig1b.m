


clc
close all


w_es = readmatrix('fig1b.txt');
ww = real(w_es(:,1));
es_iw = (w_es(:,2:8));

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

