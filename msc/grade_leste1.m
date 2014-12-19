clear all;close all;clc

load posicoes_leste1.dat;

pb=posicoes_leste1;

nest=pb(:,1);
lat=pb(:,2);
lon=pb(:,3);

lonlim=[-42 -32];
latlim=[-21 -9];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

[data,lon_bat,lat_bat]=m_tbase([lonlim(1) lonlim(2) latlim(1) latlim(2)]);
data2=smoo2(data,-9999,9,3.5);
load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap;

bat1=[-100:-50:-6000];
bat2=[-200 -1000 -5000];

figure(1)
  set(gcf,'Color',[1 1 1]) 
  [c,h]=m_contourf(lon_bat,lat_bat,data2,bat1); hold on
  caxis([bat1(end) 2000]); shading flat; 
  colormap(ocean_colormap);
  [c,h]=m_contour(lon_bat,lat_bat,data2,bat2,'k'); 
  clabel(c,h,'VerticalAlignment','middle','fontsize',10,'labelspacing',400);
  m_usercoast('../../common/costa_leste.mat','patch',[.7 .7 .7]);hold on
k=1;
%  for i=[1:51 53:63 65 67 68 70:80 82:113]
%     est=num2str(i);
%     eval(['m_text(lon(k),lat(k),est)']);
%     k=k+1;
%  end
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
  m_plot(lon,lat,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',4,'MarkerFaceColor',...
   'r','MarkerEdgeColor','w');
  m_text(-41,-10,'Grade Leste I','color','y','fontsize',10,'fontweight','bold');

 print -depsc ../figuras/grade_leste1.eps
!epstopdf ../figuras/grade_leste1.eps
%  !rm -rf ../figuras/grade_leste1.eps
