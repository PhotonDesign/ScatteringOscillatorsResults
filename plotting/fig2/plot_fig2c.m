
% *** plot fig2c
clc
close all


w_htc = readmatrix('fig2c.txt');
freq = w_htc(:,1);
bound_whtc = w_htc(:,2);
optimal_whtc = w_htc(:,3);

% *** plot setup
wToeV = 6.5821e-16;
cbound = [0.7,0,0]; coptimal = [0.3477, 0.3984, 0.7188];

% *** plot spectral HTC
lw = 3; fs = 37;
figure1 = figure;
set(figure1,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
set(gcf,'color','w');
axes1 = axes('Parent',figure1);
set(axes1,'fontsize',fs)
hold(axes1,'on');
box(axes1,'on');

p2=semilogy(freq,bound_whtc,'linewidth',lw,'color',cbound);
p1=semilogy(freq,optimal_whtc,'linewidth',lw,'color',coptimal);
hold on
set(gca, 'YScale', 'log')
xlim([0.01,0.2])
set(gca,'xaxislocation','origin','fontsize',fs-5)
ax = gca; ax.XAxis.Exponent = 0;
hold on 
box off
ylabel('Spectral HTC (Wm^{-2}K^{-1}eV^{-1})','interpreter','tex','fontsize',fs+3);
xlabel('Frequency, \omega (eV)','fontsize',fs+3,'interpreter','tex'); 
set(gcf,'position',[10,200,900,800]); 

annotation('textbox',[0.3,0.5,0.3,0.12],'EdgeColor','none','String','ideal metal','color',coptimal,'FontSize',fs-5,'Linewidth',lw-1)
annotation('textbox',[0.45,0.67,0.3,0.12],'EdgeColor','none','String','bound','color',cbound,'FontSize',fs-5,'Linewidth',lw-1)

set(gcf,'renderer','painters')
print -depsc2 spectral.eps