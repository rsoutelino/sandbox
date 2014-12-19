%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%     ANALISE DE MASSAS DE AGUA DA OPERACAO OCEANO LESTE I      %
%                 (OE-I) DA MARINHA DO BRASIL                   %
%                                                               %
%   Rafael Augusto de Mattos                                    %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%   Determinacao dos indices termohalinos e interfaces
%   entre as massas de agua ACAS, AIA e APAN
%

clc;close all;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Leitura dos Perfis Medios de Salinidade, %%%
%%%     Temperatura, Pressao e Densidade     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load medios;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Diagram T-S Medio %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[Sg,Tg] = meshgrid(32:.25:38,-2:30);
Dg = sw_dens0(Sg,Tg)-1000;

figure(1)
[c,h] = contour(Sg,Tg,Dg,20:1:40,'k');
set(h,'Color',[0.7 0.7 0.7]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
hold on;
plot(Sm,Tpm,'b.','MarkerSize',2);
title('Triangulo de Mistura','fontweight','bold','fontsize',12)
xlabel('Salinidade','fontsize',12)
ylabel('Temperatura Potencial [ \circ C]','fontsize',12)
drawnow
set(1,'Color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indice Termohalino ACAS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  acas = [36.41 20];   % <-- Silva,1995
%  acas = [36.36 20]   % <-- Miranda,1985
%  acas = [36.48 20.72]   % <-- BIO
%  acas = [36.47 20.85]   % <-- BIO I
%  acas = [36.52 20.88]   % <-- BIO II

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inicio do Processo de Analise %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot1 = plot(acas(1),acas(2),'ko','MarkerFaceColor','k','MarkerSize',7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinacao da Reta de Mistura (ACAS-AIA) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perg = 1;

while perg == 1
  aia = ginput(1);
  %%% Ajuste de Reta --> T(S) = AS + B %%%
  coef1 = polyfit([aia(1) acas(1)],[aia(2) acas(2)],1);
  plot2 = plot(aia(1),aia(2),'ko','MarkerFaceColor','k','MarkerSize',7);
  plot3 = plot([acas(1) aia(1)],[acas(2) aia(2)],'k-','linewidth',2);
  perg = menu('Manter reta de mistura ACAS-AIA ?','NAO','SIM');
    if perg == 1
      set(plot2,'Marker','none');
      set(plot3,'LineStyle','none')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinacao da Reta de Mistura (AIA-APAN) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perg = 1;

while perg == 1
  apan1 = ginput(1);
  plot4 = plot(apan1(1),apan1(2),'ko','MarkerFaceColor','k','MarkerSize',7);
  apan2 = ginput(1);
  plot5 = plot(apan2(1),apan2(2),'ko','MarkerFaceColor','k','MarkerSize',7);
  plot6 = plot([apan1(1) apan2(1)],[apan1(2) apan2(2)],'k-','linewidth',2);
  %%% Ajuste de Reta --> T(S) = AS + B %%%
  coef2 = polyfit([apan1(1) apan2(1)],[apan1(2) apan2(2)],1);
  perg = menu('Manter reta de mistura AIA-APAN ?','NAO','SIM');
  if perg == 1
      set(plot4,'Marker','none');
      set(plot5,'Marker','none');
      set(plot6,'LineStyle','none')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinacao da Reta de Mistura (APAN-AAF) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perg = 1;

while perg == 1
  aaf1 = ginput(1);
  plot7 = plot(aaf1(1),aaf1(2),'ko','MarkerFaceColor','k','MarkerSize',7);
  aaf2 = ginput(1);
  plot8 = plot(aaf2(1),aaf2(2),'ko','MarkerFaceColor','k','MarkerSize',7);
  plot9 = plot([aaf1(1) aaf2(1)],[aaf1(2) aaf2(2)],'k-','linewidth',2);
  %%% Ajuste de Reta --> T(S) = AS + B %%%
  coef3 = polyfit([aaf1(1) aaf2(1)],[aaf1(2) aaf2(2)],1);
  perg = menu('Manter reta de mistura APAN-AAF ?','NAO','SIM');
  if perg == 1
      set(plot7,'Marker','none');
      set(plot8,'Marker','none');
      set(plot9,'LineStyle','none')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinando pontos de interseccao entre as retas (Indices Termohalinos) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indice Termohalino AIA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Saia = (coef2(2)-coef1(2))/(coef1(1)-coef2(1));
Taia = coef1(1)*Saia+coef1(2);
aia = [Saia Taia];   % <-- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indice Termohalino APAN %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sapan = (coef3(2)-coef2(2))/(coef2(1)-coef3(1));
Tapan = coef2(1)*Sapan+coef2(2);
apan = [Sapan Tapan];   % <-- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotagem do Triangulo de Mistura ACAS-AIA-APAN %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(plot1,'Marker','none');
set(plot2,'Marker','none');
set(plot3,'Linestyle','none');
set(plot4,'Marker','none');
set(plot5,'Marker','none');
set(plot6,'Linestyle','none');
set(plot7,'Marker','none');
set(plot8,'Marker','none');
set(plot9,'LineStyle','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lados do Triangulo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot10 = plot([acas(1) aia(1) apan(1) acas(1)],[acas(2) aia(2) apan(2) acas(2)],...
              'k-*','linewidth',1,'MarkerEdgeColor','r');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mediana Principal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

plot11 = plot([aia(1) mean([acas(1) apan(1)])],[aia(2) mean([acas(2) apan(2)])],...
              'k--','linewidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Medianas Secundarias %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot12 = plot([mean([aia(1) acas(1)]) mean([acas(1) apan(1)])],...
              [mean([aia(2) acas(2)]) mean([acas(2) apan(2)])],...
              'r--','linewidth',1);
plot13 = plot([mean([aia(1) apan(1)]) mean([acas(1) apan(1)])],...
              [mean([aia(2) apan(2)]) mean([acas(2) apan(2)])],...
              'r--','linewidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determinacao das interfaces %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Interface AT-ACAS %%%

dens0 = sw_dens0(acas(1),acas(2))-1000;
pint0 = pm(near(Dm,dens0,1));

%%% Interface ACAS-AIA %%%

afr = 1;

menu(['Clique em OK para selecionar o ponto de interseccao'...
     ' entre a primeira mediana secundaria e a curva T-S'],'OK');
while afr == 1
  inter1 = ginput(1);
  dens1 = sw_dens0(inter1(1),inter1(2))-1000;

  ff = find(Dm == dens1);
  if isempty(ff) == 1
    afr = menu(['O ponto selecionado nao se encontra sobre a curva T-S.'...
               ' Deseja selecionar novamente?'],'SIM','NAO');
    if afr == 2
      diffd = abs(Dm-dens1);
      pint1 = pm(find(diffd == min(diffd)));
    end
  elseif isempty(ff) == 0
    pint1 = pm(ff);
    afr = 2;
  end
end

%%% Interface AIA-APAN %%%

afr = 1;

menu(['Clique em OK para selecionar o ponto de interseccao'...
     ' entre a segunda mediaana secundaria e a curva T-S'],'OK');
while afr == 1
  inter2 = ginput(1);
  dens2 = sw_dens0(inter2(1),inter2(2))-1000;
  
  ff = find(Dm == dens2);
  if isempty(ff) == 1
    afr = menu(['O ponto selecionado nao se encontra sobre a curva T-S.'...
               ' Deseja selecionar novamente?'],'SIM','NAO');
    if afr == 2
      diffd = abs(Dm-dens2);
      pint2 = pm(find(diffd == min(diffd)));
    end
  elseif isempty(ff) == 0
    pint2 = pm(ff);
    afr = 2;
  end
end

%%% ISOPICNAIS DAS INTERFACES NO DIAGRAMA T-S %%%

[c,h] = contour(Sg,Tg,Dg,[dens0 dens1 dens2],'m');
cl = clabel(c,h,'manual','VerticalAlignment','middle','fontsize',10,'labelspacing',500);
for j = 1:length(h)
  strh = get(cl(j),'String');strh = strh(1:5);
  set(cl(j),'String',strh,'Color','m');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indices termohalinos %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Tacas = num2str(acas(1,2),4); Sacas = num2str(acas(1,1),4);
Taia = num2str(aia(1,2),3); Saia = num2str(aia(1,1),4);
Tapan= num2str(apan(1,2),3); Sapan = num2str(apan(1,1),4);

text(acas(1,1)-0.2,acas(1,2),['ACAS (',Sacas,',',Tacas,')'],'BackgroundColor',[1 1 1],...
     'Margin',3,'Fontweight','bold','HorizontalAlignment','right');
text(aia(1,1)-0.01,aia(1,2)-1,['AIA (',Saia,',',Taia,')'],'BackgroundColor',[1 1 1],...
     'Margin',3,'Fontweight','bold','HorizontalAlignment','right');
text(apan(1,1)+0.2,apan(1,2),['APAN (',Sapan,',',Tapan,')'],'BackgroundColor',[1 1 1],...
     'Margin',3,'Fontweight','bold','HorizontalAlignment','left');

hold off;

%  saveas(1,'figuras_rad2/triangulo_rad2','fig');
%  print(1,'-depsc','figuras_rad2/triangulo_rad2')
%  eval(['!epstopdf figuras_rad2/triangulo_rad2.eps']);
%  eval(['!rm figuras_rad2/triangulo_rad2.eps']);

disp(' ');
disp([' Interface AT-ACAS = [',num2str(pint0),' m; ',...
     num2str(dens0),' kgm-3; ',num2str(Tm(find(pm == pint0))),' C; ',num2str(Sm(find(pm == pint0))),']']);
disp(' ');
disp([' Interface ACAS-AIA = [',num2str(pint1),' m; ',...
     num2str(dens1),' kgm-3; ',num2str(Tm(find(pm == pint1))),' C; ',num2str(Sm(find(pm == pint1))),']']);
disp(' ');
disp([' Interface AIA-APAN = [',num2str(pint2),' m; ',...
     num2str(dens2),' kgm-3; ',num2str(Tm(find(pm == pint2))),' C; ',num2str(Sm(find(pm == pint2))),']']);
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interfaces no Perfil de SigmaTheta %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(Dm,-pm,'m','linewidth',1);
hold on;

plot14 = plot([get(gca,'XLim')],-[pint0 pint0],'r--','linewidth',2);

plot15 = plot([get(gca,'XLim')],-[pint1 pint1],'r--','linewidth',2);

plot16 = plot([get(gca,'XLim')],-[pint2 pint2],'r--','linewidth',2);

text(25.1,-pint0-90,['interface AT-ACAS = -',num2str(pint0),' m'],'BackgroundColor',[1 1 1]);
text(25.1,-pint1-90,['interface ACAS-AIA = -',num2str(pint1),' m'],'BackgroundColor',[1 1 1]);
text(25.1,-pint2-90,['interface AIA-APAN = -',num2str(pint2),' m'],'BackgroundColor',[1 1 1]);
hold off;
set(gca,'YLim',[-1400 0]);
title('Perfil vertical medio de \sigma_{\theta}','fontsize',12,'fontweight','bold')
xlabel('Densidade Potencial [kg m^{-3}]','fontsize',12)
ylabel('Profundidade [m]','fontsize',12)
drawnow;
set(2,'Color','w')

 saveas(1,'triangulo','fig');
print(1,'-depsc','triangulo');
!epstopdf triangulo.eps
!rm -rf triangulo.eps
print(2,'-depsc','dens_interf');
!epstopdf dens_interf.eps
!rm -rf dens_interf.eps


