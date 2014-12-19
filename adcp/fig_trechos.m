clear all;close all;clc;warning off

load contour/oeii_full_xy.mat;

lon = xyt(1,:); lon = lon-360;
lat = xyt(2,:);
t = xyt(3,:);

lonlim=[-44 -32];
latlim=[-24 -9];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');
[x,y] = m_ll2xy(lon,lat);

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
  m_usercoast('costa_adcp.mat','patch',[.7 .7 .7]);hold on
k=1;
%  for i=[1:4 7:112]
%     est=num2str(i);
%     eval(['m_text(lon(k),lat(k),est)']);
%     k=k+1;
%  end
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%    m_plot(lon,lat,'Color','r','LineStyle','none','LineWidth',...
%         1,'marker','o','markersize',4,'MarkerFaceColor',...
%     'r','MarkerEdgeColor','w');
  m_text(-43,-10,'Trajeto do ADCP - OEII','color','y','fontsize',10,'fontweight','bold');

arrow([x(1) y(1)],[x(11) y(11)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(11) y(11)],[x(21) y(21)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(21) y(21)],[x(31) y(31)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(31) y(31)],[x(80) y(80)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(81) y(81)],[x(121) y(121)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(122) y(122)],[x(158) y(158)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(158) y(158)],[x(174) y(174)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(174) y(174)],[x(213) y(213)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(216) y(216)],[x(250) y(250)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(250) y(250)],[x(256) y(256)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(256) y(256)],[x(290) y(290)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(296) y(296)],[x(330) y(330)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(330) y(330)],[x(346) y(346)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(346) y(346)],[x(392) y(392)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(445) y(445)],[x(481) y(481)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(481) y(481)],[x(491) y(491)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(491) y(491)],[x(524) y(524)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(525) y(525)],[x(554) y(554)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(554) y(554)],[x(562) y(562)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(562) y(562)],[x(588) y(588)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(589) y(589)],[x(620) y(620)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(620) y(620)],[x(629) y(629)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
arrow([x(629) y(629)],[x(673) y(673)],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')
%  arrow([x() y()],[x() y()],'length',7,'tipangle',5,'linewidth',1,'facecolor','y','edgecolor','y')

 print -depsc ../figuras/adcp_trechos.eps
!epstopdf ../figuras/adcp_trechos.eps
%  !rm -rf ../figuras/grade_leste2.eps
