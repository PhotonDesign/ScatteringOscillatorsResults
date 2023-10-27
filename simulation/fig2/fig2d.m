%% *** plot fig2d
clc
close all

ms = '*'; msize = 12;

lw = 3; fs = 40; cblue = [0 0.4470 0.7410]./1;
kb = 1.3806e-23;
hbar = 1.0546e-34;
q = 1.6022e-19;
T = 300; Tscaled = T * kb / hbar;
wToeV = hbar/q;
d = 1e-8;

c2D = 0.9.*[0.0195, 0.0859, 0.75]; 
coptimal_bulk = 0.8.*[0.9961, 0.0078, 0.0078];
cbound = [1,1,1].*0.4;

TT = 100.*linspace(1,12,12);
T_compute = [100, 300, 600, 900, 1200];
wien_compute = [0.02272, 0.0683, 0.137, 0.2072, 0.273];
wien_optimal = interp1(T_compute,wien_compute,TT,'linear')./wToeV;
TT_s = TT.*kb./hbar;
w_0_nf_htc = 2.57368.*TT_s;

figwien = figure;
set(gcf,'color','w');
set(figwien,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
axes1 = axes('Parent',figwien);
set(axes1,'fontsize',fs')

%yyaxis left
plot(TT,w_0_nf_htc,'linewidth',lw,'color',cbound)
hold on
plot(TT,wien_optimal,'linestyle','none','linewidth',lw-1.5,'markersize',msize,'Marker',ms,'color',c2D,'markerfacecolor',c2D)
plot(TT,w_0_nf_htc,'linestyle','none','linewidth',lw-1.5,'markersize',msize,'Marker',ms,'color',coptimal_bulk,'markerfacecolor',coptimal_bulk)
xlabel('Temperature, T (K)','fontsize',fs+3,'interpreter','tex');
ylabel('Optimal frequency, \omega (eV)','fontsize',fs+3,'interpreter','tex','position',[-120,1.9e14]); 
xlim([0,1200])
xticks([0,300,600,900,1200])
wtick = linspace(0.00,0.3,4)./wToeV;
yticks(wtick);
yticklabels({'0','0.1','0.2','0.3','0.4'});
ylim([0,0.3]./wToeV)
hold on
set(gca,'fontsize',fs) 
set(gcf,'position',[10,200,800,600]); 
box off

set(gcf,'renderer','painters')
print -depsc2 wien_bulk_2D.eps