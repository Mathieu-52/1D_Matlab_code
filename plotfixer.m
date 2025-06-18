xlabel('$x$','interpreter','latex','fontsize',16);%'fontsize',48,'Fontweight','bold')
%ylabel('$$','interpreter','latex','fontsize',15);%'fontsize',48,'Fontweight','bold')
%ylabel('Final Error','fontsize',48,'Fontweight','bold')
axis tight;
%plot([-0.04 0.12],[0.02 -0.06]-0.015,'--k','linewidth',2)

set(gca,'Fontsize',16);
% h_legend=legend('${h}^{+}/R,U=0.5m/s,R=0.85mm$','${h}^{-}/R,U=0.5m/s,R=0.85mm$','${h}^{+}/R,U=0.4m/s,R=1.15mm$','${h}^{-}/R,U=0.4m/s,R=1.15mm$','${h}^{+}/R,U=0.7m/s,R=0.55mm$','${h}^{-}/R,U=0.7m/s,R=0.55mm$','$-1/2\: slope$');
% set(h_legend,'interpreter','latex','Fontsize',16);
% legend('boxoff')

file_name=sprintf('Paper_Figures/eqbm_D1_1_D2_1_T_timelapse');

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
saveas(gcf,file_name,'fig') 
print(gcf,'-depsc','-r100',file_name);