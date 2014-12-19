clear all; close all; clc;warning off

load ../proc/hidrografia/posicoes_leste2.dat;
pb=posicoes_leste2;
nest=pb(:,1);
latg=pb(:,2);
latm=pb(:,3);
long=pb(:,4);
lonm=pb(:,5);
prof=pb(:,6);
lat=-(latg+latm/60);
lon=-(long+lonm/60);
lonlim=[-42 -30];
latlim=[-21 -9];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

[data,lon_bat,lat_bat]=m_tbase([lonlim(1) lonlim(2) latlim(1) latlim(2)]);
data2=smoo2(data,-9999,9,3.5);
bat2=[-200 -1000 -5000];

figure(1)
  set(gcf,'Color',[1 1 1]) 
  [c,h]=m_contour(lon_bat,lat_bat,data2,bat2,'k'); clabel(c,h,'VerticalAlignment','middle','fontsize',10,'labelspacing',400);
  m_usercoast('../proc/common/costa_leste.mat','patch',[0 0 0]);hold on
 m_plot(lon,lat,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',4,'MarkerFaceColor',...
   'r','MarkerEdgeColor','w');
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',14);

sup = [0.0411 -0.2038;
       0.0209 -0.2368;
      -0.0342 -0.2183;
      -0.0116 -0.1835]; xsup = sup(:,1); ysup = sup(:,2);
cen = [0.0162 -0.2386;
       0.0156 -0.3023;
      -0.0556 -0.3023;
      -0.0539 -0.2357]; xcen = cen(:,1); ycen = cen(:,2);
inf = [0.0220 -0.3127;
       0.0220 -0.3521;
      -0.0666 -0.3521;
      -0.0620 -0.3116]; xinf = inf(:,1); yinf = inf(:,2);
fill(xsup,ysup,[1 0.951 0.495]);
fill(xcen,ycen,[0.843 1 0.335]);
fill(xinf,yinf,[0.614 1 1]);

print -depsc ../figuras/setores.eps
























