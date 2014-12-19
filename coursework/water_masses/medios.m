%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Lista - Sueli
%%%    Oceanografia Regional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

Tm = [];Sm = [];pm = []; lon=[]; lat=[];

xx = 1:12; 
Se = 33.5:.25:38;
Te = -2:30;

[Sg,Tg] = meshgrid(Se,Te);
dens = sw_dens0(Sg,Tg)-1000;

%  figure(1)
%  [c,h] = contour(Se,Te,dens,20:1:40,'k');
%  cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
%  xlabel('Salinidade','fontsize',12)
%  ylabel('Temperatura [ \circ C]','fontsize',12)
%  hold on
%  drawnow
%  set(1,'Color','w')

figure(10)
[c,h] = contour(Se,Te,dens,20:1:40,'k');
set(h,'Color',[0.7 0.7 0.7]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
title('Diagrama de Estado','fontweight','bold','fontsize',12)
xlabel('Salinidade','fontsize',12)
ylabel('Temperatura [ \circ C]','fontsize',12)
hold on
drawnow
set(10,'Color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOPING PRINCIPAL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:length(xx)
  
  XX = num2str(xx(k));
  eval(['load dados/est',XX,'.dat;']);
  eval(['ctd = est',XX,';']);
  eval(['clear est',XX]);
  p = ctd(2:end,1);
  T = ctd(2:end,2);
  S = ctd(2:end,3);
  lon = [lon; ctd(1,3)];
  lat = [lat; ctd(1,2)];
  
  T = T(find(p>=15));
  S = S(find(p>=15));
  p = p(find(p>=15));

  
    if size(T,1) <= size(Tm,1)
      Tm(:,k) = NaN; Sm(:,k) = NaN; pm(:,k) = NaN;
      Tm(1:size(T,1),k) = T; Sm(1:size(S,1),k) = S; pm(1:size(p,1),k) = p;
    else
      TM = NaN*ones(size(T,1),k); SM = NaN*ones(size(S,1),k); pM = NaN*ones(size(p,1),k);
      TM(:,k) = T; SM(:,k) = S; pM(:,k) = p;
      TM(1:size(Tm,1),1:size(Tm,2)) = Tm; SM(1:size(Sm,1),1:size(Sm,2)) = Sm; pM(1:size(pm,1),1:size(pm,2)) = pm;
      Tm = TM; Sm = SM; pm = pM;
      clear TM SM pM;
    end
   
    figure(10)
    plot(S,T,'b.','markersize',0.1)

%      figure(2)
%      plot(T,-p,'b:','linewidth',.5);
%      hold on;
%      xlabel('Temperatura [ \circ C]','fontsize',12)
%      ylabel('Profundidade [m]','fontsize',12)
%      set(2,'Color','w')
%    
%      figure(3)
%      plot(S,-p,'b:','linewidth',.5);
%      hold on;
%      xlabel('Salinidade','fontsize',12)
%      ylabel('Profundidade [m]','fontsize',12)
%      set(3,'Color','w')
  end
    
  clear p T S XX ctd
   


%%% PERFIS MEDIOS %%%

Tm(find(Tm == 0)) = NaN;
Sm(find(Sm == 0)) = NaN;
pm(find(pm == 0)) = NaN;

pm = nanmean(pm')';
Tm = nanmean(Tm')';
Sm = nanmean(Sm')';
Dm = sw_pden(Sm,Tm,pm,0)-1000;
Tpm = sw_ptmp(Sm,Tm,pm,0);

%  save medios pm Tm Sm Dm Tpm

atblocX = [36 36 37.5 37.5 36];
atblocY = [28 20.2 20.2 28 28];

acasblocX = [35 35 36.4 36.4 35];
acasblocY = [20 10.84 10.84 20 20];

aiablocX = [34 34 35.2 35.2 34];
aiablocY = [10.8 4 4 10.8 10.8];

apanblocX = [34.2 34.2 35 35 34.2];
apanblocY = [3.83 2 2 3.83 3.83];

figure(10)
plot(atblocX,atblocY,'m','linewidth',1.5)
gtext('AT','color','m')
plot(acasblocX,acasblocY,'g','linewidth',1.5)
gtext('ACAS','color','g')
plot(aiablocX,aiablocY,'k','linewidth',1.5)
gtext('AIA','color','k')
plot(apanblocX,apanblocY,'r','linewidth',1.5)
gtext('APAN','color','r')

print(10,'-depsc',['TS_espalhado']);
!epstopdf TS_espalhado.eps
!rm -rf TS_espalhado.eps


