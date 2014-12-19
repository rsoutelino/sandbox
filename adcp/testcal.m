clear all;close all;clc

cruz = 'exper1';
CRUZ = 'OEII';

eval(['load ../',cruz,'/contour/contour_uv.mat'])
eval(['load ../',cruz,'/contour/contour_xy.mat'])
u = uv(:,1:2:end-1);
v = uv(:,2:2:end);
lon = xyt(1,:)-360; lat = xyt(2,:);
t = xyt(3,:);
eval(['load ../',cruz,'/cal/watertrk/oeii_7.cal']);
eval(['tcal = oeii_7(:,2);']);
eval(['ang = oeii_7(:,11);']); 

F = [];

for k = 1:length(tcal)
    f = near(t,tcal(k),1);
    F = [F f];
end

loncal = lon(F); latcal = lat(F);

%% pegando o angulo medio de fase PHI
eval(['fid = fopen(''../',cruz,'/cal/watertrk/adcpcal.out'',''r'') '])
for i = 1:12; fgetl(fid); end
linha = fgetl(fid);
phi = linha(22:25);

%% configurando m_map

load ../../common/etopo2_leste.mat

figure;
set(gcf,'color','w');hold on
  lonlim=[-41 -33.7]; latlim=[-20.5 -10.2];
  m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');
  bat1 = [0:100:5000];
  bat2 = [200 1000 2500 3000];
  load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap
  cm1 = colormap(ocean_colormap);cm1 = flipud(cm1);
  m_contourf(xb,yb,-zb,bat1); shading flat; colormap(cm1);%colorbar
  [c,h] = m_contour(xb,yb,-zb,bat2,'k');hold on
  c = clabel(c,h,'labelspacing',1000);
  set(h,'color',[.4 .4 .4])
  set(c,'color',[.4 .4 .4],'fontsize',8)

  m_usercoast('../../common/costa_leste.mat','patch',[0 0 0],'LineStyle','-');
  m_grid('box','fancy')

m_plot(lon,lat,'w');hold on
%  m_plot(lon,lat,'k.','markersize',0.9)
m_plot(loncal,latcal,'yo')
m_plot(-40.5,-11,'yo'); m_text(-40.3,-11,'Pontos de Calibracao','color','w','fontsize',8,'fontweight','bold')
t = ['\lambda Medio = ',phi,' ^{\circ}'];
m_text(-40.5,-11.3,t,'color','w','fontweight','bold','fontsize',8)
tit = ['Pontos de Calibracao - CODAS - ',CRUZ];
title(tit,'fontweight','bold')

eval(['print -depsc ../figuras/calib_codas_',CRUZ,'.eps']);
eval(['!epstopdf ../figuras/calib_codas_',CRUZ,'.eps']);






















