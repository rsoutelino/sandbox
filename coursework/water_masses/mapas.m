%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Rotina para plotar mapas         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

latlim=[-27 -23];
lonlim=[-48 -43];

%limmer=['   Limites Zonais: ',num2str(lonlim(1)),' ',num2str(lonlim(2))];
%disp(limmer);
%limzon=['   Limites Meridionais: ',num2str(latlim(1)),' ',num2str(latlim(2))];
%disp(limzon);
%disp(' ');
%disp('Obs.: Longitude Oeste [W] possui sinal NEGATIVO - Longitude Leste [E] possui sinal POSITIVO');
%disp('      Latitude Sul [S] possui sinal NEGATIVO - Latitude Norte [N] possui sinal POSITIVO');
%disp(' ');

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

[data,lon_bat,lat_bat]=m_tbase([lonlim(1) lonlim(2) latlim(1) latlim(2)]);
data2=smoo2(data,-9999,9,3.5);

disp('Determinando batimetria e linha de costa...')

[lon_bat,lat_bat]=m_ll2xy(lon_bat,lat_bat);

load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap;
colormap(ocean_colormap);

bat1=[-100:-50:-6000];
bat2=[-200 -1000 -5000];


figure(1); axes('color','none'); hold on;

[c,h]=contourf(lon_bat,lat_bat,data2,bat1);
%  caxis([bat1(end) 2000]);
caxis([-2500 0]);
shading flat;hold on;
[c,h]=contour(lon_bat,lat_bat,data2,bat2,'k'); 
clabel(c,h,'VerticalAlignment','middle','fontsize',10,'labelspacing',400);

m_gshhs_i('patch',[.6 .6 .6]);  
m_grid('box','fancy','xtick',4,'ytick',5,'yaxislocation','right','xaxislocation','top','fontsize',12);
m_grid('box','fancy','xtick',4,'ytick',5,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
grid off
set(gcf,'color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% PLOT RADIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load posicoes.mat; 

m_plot(lon,lat,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',5,'MarkerFaceColor',...
       'r','MarkerEdgeColor','w');

%m_text(lon,lat+.02,num2str(nest),'color','y','vertical',...
%       'bottom','FontWeight','bold','fontsize',10);

print(1,'-depsc','mapa.eps')
!epstopdf mapa.eps
!rm -rf  mapa.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
