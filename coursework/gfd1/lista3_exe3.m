clear all;close all;clc

x = 1:400;
y = 1:20;
k = 0.0628;

eta = cos(k*x);
ETA = (ones(size(y'))*eta);

yy = (ones(size(x'))*y)';



figure(1)
set(gcf,'color','w')
contour(x,y,ETA)
cc = colorbar
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*2.7 pos(3)/4 pos(4)/1.8])
%  plot(hlambda,c_cl,'linewidth',2); hold on
%  plot(hlambda,c_cc,'r','linewidth',2)
%  axis([hlambda(1) hlambda(end)+0.05 0 1.2])
%  plot(hlambda,0.95*ones(size(hlambda)),'--k');
%  legend('C/Cl','C/Cc','0.95',0)
title('Isolinhas de \eta [m] em t=0','fontweight','bold')
ylabel('y [m]','fontweight','bold')
xlabel('x [m]','fontweight','bold')
pbaspect([1 0.5 1])

print -depsc exe3.eps 
!epstopdf exe3.eps
!rm -rf exe3.eps

figure(2)
set(gcf,'color','w')
surf(x,y,ETA);shading flat; hold on
cc = colorbar
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*2.7 pos(3)/4 pos(4)/1.8])
%  plot(hlambda,c_cl,'linewidth',2); hold on
%  plot(hlambda,c_cc,'r','linewidth',2)
%  axis([hlambda(1) hlambda(end)+0.05 0 1.2])
%  plot(hlambda,0.95*ones(size(hlambda)),'--k');
%  legend('C/Cl','C/Cc','0.95',0)
title('Estrutura tridimensional da onda no canal','fontweight','bold')
ylabel('y [m]','fontweight','bold')
xlabel('x [m]','fontweight','bold')
zlabel('\eta [m]','fontweight','bold')
pbaspect([1 0.2 0.2])

print -depsc exe3_fig2.eps 
!epstopdf exe3_fig2.eps
!rm -rf exe3_fig2.eps


