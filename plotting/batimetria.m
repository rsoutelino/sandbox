%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PROGRAMA PARA PLOTAR A BATIMETRIA 
%	   DO CANAL DE SAO SEBASTIAO
%	   maio/2006 - Observacional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all; clc

%%% carregando os dados de posicao das estacoes
load estacoes.dat;
xest=estacoes(:,1);
yest=estacoes(:,2);

%%% carregando a batimetria e separando as colunas com lon, lat, e prof
load batimetria_GMS.dat;
ssb = batimetria_GMS;
lon=ssb(:,1);
lat=ssb(:,2);
z=ssb(:,3);

%%% estabelecendo limites para a projeçao
lonlim=[min(lon) max(lon)];
latlim=[min(lat) max(lat)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

%%% CRIANDO GRADE PARA INTERPOLAR %%%%%%%%

% grade regular para interpolar batimetria (dados estavam desordenados)
x=min(min(lon)):0.0003:max(max(lon)); % resolucao de 0.0003º
y=min(min(lat)):0.0003:max(max(lat));

[xc,yc]=meshgrid(x,y);
disp('Interpolando a Batimetria')
[x,y,zc]=griddata(lon,lat,z,xc,yc); 

%  PLOTANDO BATIMETRIA 3D

figure(1)
  set(gcf,'Color',[1 1 1]);

surf(xc,yc,zc);shading flat;
title('Batimetria - Canal de Sao Sebastiao','fontsize',12,'fontweight','bold')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Profundidade')
%  
%   print -dpng batimetria3D.png

%%% PLOTANDO BATIMETRIA NO MAPA, COM AS ESTACOES HIDROGRAFICAS

figure(2)
  set(gcf,...
        'Color',[1 1 1],... 
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

  m_contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat;
  hold on
  colorbar('horiz')
  m_usercoast('../ctd/costa_ssebastiao.mat','patch',[.6 .6 .6],'LineStyle','-');
%    m_gshhs_f('patch',[.6 .6 .6],'LineStyle','-');
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',12);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
  m_plot(xest,yest,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',8,'MarkerFaceColor',...
       'k','MarkerEdgeColor','w');

 print -depsc batimetria.eps
!epstopdf batimetria.eps



