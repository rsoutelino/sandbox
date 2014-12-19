%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% item 7 lista da Sue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

load medios.mat
load coef.mat

fs = find(Sm >= 33.2 & Sm <= 36.41);
Tm = Tm(fs);
Sm = Sm(fs);
pm = pm(fs);
a = a(fs);
b = b(fs);

ft = find(Tm >= 4 & Tm <= 20);
Tm = Tm(ft);
Sm = Sm(ft);
pm = pm(ft);
a = a(ft);
b = b(ft);

xx = 1:12; 
Se = 33.5:.25:38;
Te = -2:30;

[Sg,Tg] = meshgrid(Se,Te);
dens = sw_dens0(Sg,Tg)-1000;

figure(1)
[c,h] = contour(Se,Te,dens,20:1:40,'k');
set(h,'Color',[0.7 0.7 0.7]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
title('Curva de Mistura - ACAS/AIA','fontweight','bold','fontsize',12)
xlabel('Salinidade','fontsize',12)
ylabel('Temperatura [ \circ C]','fontsize',12)
hold on
plot(Sm,Tm,'b','markersize',1.5)
drawnow
set(1,'Color','w')

print -depsc curva_acas_aia.eps
!epstopdf curva_acas_aia.eps
! rm -rf  curva_acas_aia.eps

%%% calculando derivada termohalina

dSdT = diff(Sm)./diff(Tm);

mix = (b(2:end)./a(2:end)).*dSdT; 

figure(2)

plot([0 0],[-100 -1000],'k');hold on
plot([1 1],[-100 -1000],'b--')
plot(mix,-pm(2:end),'r','linewidth',1.5);
set(2,'Color','w')
title('Derivada Temohalina ao longo da curva','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
print -depsc dSdT.eps
!epstopdf dSdT.eps
! rm -rf  dSdT.eps




