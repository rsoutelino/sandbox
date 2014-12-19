% plota perfis hidrograficos para visualizacao rapida
% Oceano Leste I - Mestrado

clear all;close all;clc
T=[]; S=[]; 
cont = 0;
for i=[1:51 53:63 65 67 68 70:80 82:113]
      eval(['load ../../../dados/leste1/filtrados/lesteI_ctd',num2str(i),'.dat;'])
      eval(['ctd = lesteI_ctd',num2str(i),';'])
      eval(['clear lesteI_ctd',num2str(i)])

      p = ctd(:,1);
      t = ctd(:,2);
      s = ctd(:,3);
      T = [T;t];
      S = [S;s];
      
      if max(p) > 2400
         f = find(p>=15 & p<=2399);
         tm = t(f);
         sm = s(f);
         if i == 1
            Tm = tm;
            Sm = sm;
            cont = cont+1;
         else
            Tm = Tm + tm;
            Sm = Sm + sm;
            cont = cont+1;
         end
      end   
      
%  figure(1)
%     subplot(121)
%        plot(t,-p)
%        tit=['Estacao',num2str(i)];
%        title(tit)
%     subplot(122)
%        plot(s,-p,'r')
%  pause
end

Tm = Tm./cont;
Sm = Sm./cont;

[Sg,Tg] = meshgrid(33.5:.25:38,0:30);
Dg = sw_dens0(Sg,Tg)-1000;

figure(1)
set(gcf,'color','w')
[c,h] = contour(Sg,Tg,Dg,20:1:40,'k');
set(h,'Color',[0.85 0.85 0.85]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
hold on;
p = plot(S,T,'k.','MarkerSize',0.1);
set(p,'Color',[0.7 0.7 0.7]);
plot(Sm,Tm,'k','linewidth',1.5)
title('T-S Diagram','fontweight','bold','fontsize',12)
xlabel('Salinity','fontsize',12)
ylabel('Potencial Temperature [ \circ C]','fontsize',12)

print -depsc TS_leste1.eps
!epstopdf TS_leste1.eps
!rm -rf TS_leste1.eps
