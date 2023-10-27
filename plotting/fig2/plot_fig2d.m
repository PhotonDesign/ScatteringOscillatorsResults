%% *** plot fig2d
clc
close all

T_wien = readmatrix('fig2d.txt');
TT = T_wien(:,1);
w_0_nf_htc = T_wien(:,2);
wien_optimal = T_wien(:,3);
w_0_nf_htc = T_wien(:,4);

ms = '*'; msize = 12;

lw = 3; fs = 40; cblue = [0 0.4470 0.7410]./1;

c2D = 0.9.*[0.0195, 0.0859, 0.75]; 
coptimal_bulk = 0.8.*[0.9961, 0.0078, 0.0078];
cbound = [1,1,1].*0.4;

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