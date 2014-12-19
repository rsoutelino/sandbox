%%
%%

clear all;close all;

load medios.mat

kt = 3.65e-3;
ks = 2.45e-3;
h = 236;
t = 1.5*1e7;
Z = -2000:50:2000;

z1t = (Z-h)/(sqrt(2*kt*t)); 
z2t = ((Z+h)/(sqrt(2*kt*t)));
z1s = ((Z-h)/(sqrt(2*ks*t))); 
z2s = ((Z+h)/(sqrt(2*ks*t))); 

phi1t = normcdf(z1t); 
phi2t = normcdf(z2t); 
phi1s = normcdf(z1s);
phi2s = normcdf(z2s);

T1 = 20;
T2 = 2.05;
T3 = 3.55;

S1 = 36.41;
S2 = 33.17;
S3 = 34.98;



T = 0.5*(T1+T3) + (T1-T2).*phi1t + (T2-T3).*phi2t;
S = 0.5*(S1+S3) + (S1-S2).*phi1s + (S2-S3).*phi2s;


xx = 1:12; 
Se = 32.5:.25:38;
Te = -2:30;

[Sg,Tg] = meshgrid(Se,Te);
dens = sw_dens0(Sg,Tg)-1000;

figure(1)
plot(S-.7,T-8,'r','linewidth',1.5);hold on % -.7 -8
plot(Sm,Tm,'b')
legend('Teorico','Observado',0)
[c,h] = contour(Se,Te,dens,20:1:40,'k');hold on
set(h,'Color',[0.7 0.7 0.7]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
title('Diagrama de Estado Teorico','fontweight','bold','fontsize',12)
xlabel('Salinidade','fontsize',12)
ylabel('Temperatura [ \circ C]','fontsize',12)
plot([S1 S2 S3 S1],[T1 T2 T3 T1],'k--')
plot(S1,T1,'*')
plot(S2,T2,'*')
plot(S3,T3,'*')
axis([33 37 0 23])
set(gcf,'color','w')

print -depsc item8.eps
!epstopdf item8.eps
!rm -rf item8.eps












