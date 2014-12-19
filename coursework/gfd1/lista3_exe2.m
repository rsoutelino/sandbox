clear all;close all;clc

hlambda = 0.01:0.001:0.6;

c_cl = sqrt(  ( tanh(2.*pi.*hlambda)  ) ./  ( 2.*pi.*hlambda   )  );

c_cc = sqrt ( tanh(2.*pi.*hlambda) );

figure(1)
set(gcf,'color','w')

plot(hlambda,c_cl,'linewidth',2); hold on
plot(hlambda,c_cc,'r','linewidth',2)
axis([hlambda(1) hlambda(end)+0.05 0 1.2])
plot(hlambda,0.95*ones(size(hlambda)),'--k');
legend('C/Cl','C/Cc','0.95',0)
title('Comparacao entre C/Cl e C/Cc','fontweight','bold')
xlabel('H/\lambda','fontweight','bold')
pbaspect([1 0.5 1])

print -depsc exe2.eps 
!epstopdf exe2.eps
!rm -rf exe2.eps


