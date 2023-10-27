clc
close all

w_eig = readmatrix('fig1d.txt');
ww = real(w_eig(:,1));
ImTeig = w_eig(:,2:end);

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
ng = 84;
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