%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Lista - Sueli
%%%    Oceanografia Regional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

Tm = [];Sm = [];pm = []; lon=[]; lat=[];

% usar xx=6 para o item 2 da lista

  load dados/est7.dat;
  ctd = est7;
  
  p = ctd(2:end,1);
  T_pot = ctd(2:end,2);
  S = ctd(2:end,3);
   
 T = sw_temp(S,T_pot,p,1); 
 dens = sw_dens(S,T,p);
 pden = sw_pden(S,T,p,1); 
 a = sw_alpha(S*0,T,p*0,'temp');
 b = sw_beta(S,T*0,p*0,'temp');
 adtg = sw_adtg(S,T,p);
 svel = sw_svel(S,T,p);
 k = sw_seck(S,T,p); k = 1./k;
 Cp = sw_cp(S,T,p*0); 

save coef.mat dens pden a b adtg svel k Cp

figure(2)

plot(T_pot,-p,'b','linewidth',1.5); hold on
plot(T,-p,'r','linewidth',1.5)
title('Temperatura [^\circ C]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
legend('Potencial','In Situ',0)
plot([4 4],[-2500 0],'k--'); hold on
set(2,'Color','w')
print(2,'-depsc',['temp.eps']);
!epstopdf temp.eps
!rm -rf  temp.eps


figure(3)
plot(S,-p,'r','linewidth',1);
hold on
title('Salinidade','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square
set(3,'Color','w')
print(3,'-depsc',['sal.eps']);
!epstopdf sal.eps
!rm -rf  sal.eps

figure(5)
plot(pden-1000,-p,'b','linewidth',1); hold on
plot(dens-1000,-p,'r','linewidth',1)
title('Densidade [kg m^{-3}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
legend('\sigma_{\theta}','\sigma_T',0)
set(5,'Color','w')
print(5,'-depsc',['dens.eps']);
!epstopdf dens.eps
!rm -rf  dens.eps

figure(6)
plot([0 0],[-2500 0],'k--'); hold on
plot(a*1e4,-p,'b','linewidth',1);
hold on
title('Coeficiente de Expansao Termica [^\circ C^{-1} . 10^{-4}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(6,'Color','w')
print(6,'-depsc',['alfa.eps']);
!epstopdf alfa.eps
!rm -rf  alfa.eps

figure(7)
plot(b*1e4,-p,'b','linewidth',1);
hold on
title('Coeficiente de Contracao Salina [^{-1} . 10^{-4}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(7,'Color','w')
print(7,'-depsc',['beta.eps']);
!epstopdf beta.eps
!rm -rf  beta.eps

figure(8)
plot(adtg,-p,'b','linewidth',1);
hold on
title('Gradiente Adiabatico de Temperatura [^\circ C dbar^{-1}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(8,'Color','w')
print(8,'-depsc',['adtg.eps']);
!epstopdf adtg.eps
!rm -rf  adtg.eps

figure(9)
plot(svel,-p,'b','linewidth',1);
hold on
title('Velocidade de propagacao do som [m s^{-1}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(9,'Color','w')
print(9,'-depsc',['svel.eps']);
!epstopdf svel.eps
!rm -rf svel.eps


figure(10)
plot(k*10e5,-p,'b','linewidth',1);
hold on
title('Coeficiente de Compressibilidade Barica [dbar^{-1} . 10^{-5}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(10,'Color','w')
print(10,'-depsc',['cbarica.eps']);
!epstopdf cbarica.eps
!rm -rf cbarica.eps


figure(11)
plot(Cp,-p,'b','linewidth',1);
hold on
title('Calor Especifico a pressao e volume constantes [J kg^{-1} ^\circ C^{-1}]','fontweight','bold','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
axis square; drawnow
set(11,'Color','w')
print(11,'-depsc',['cp.eps']);
!epstopdf cp.eps
!rm -rf cp.eps
