%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting HYCOM GLOBAL outputs
% Rafael Soutelino - May 2010
% rsoutelino@gmail.com
% creates average fields in feb and march 2005 for comparissom with OEII cruise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYCOM + NCODA Global 1/12 Analysis (expt_60.5) (assimilates data) Nov-2003 to Dec-2006
url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5'; 
quiv = 'y';
PRINT = 'y';
depth = 100; DEPTH = num2str(depth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% loading previously saved data in plot_hycom.m ======================================
eval(['load hycom_glb08_exp60.5_feb2005_',DEPTH,'m.mat']);
ufeb = ucum; vfeb = vcum; tfeb = tcum;

eval(['load hycom_glb08_exp60.5_mar2005_',DEPTH,'m.mat']);
umar = ucum; vmar = vcum; tmar = tcum;

% computing average fields ===========================================================
umar = squeeze(mean(umar)); vmar = squeeze(mean(vmar)); tmar = squeeze(mean(tmar));
ufeb = squeeze(mean(ufeb)); vfeb = squeeze(mean(vfeb)); tfeb = squeeze(mean(tfeb));
uoe2 = (umar + ufeb)./2; voe2 = (vmar + vfeb)./2; toe2 = (tmar + tfeb)./2;

% plotting average fields ============================================================
%  figure
%  set(gcf,'Position',[0 400 1200 650],'color','w')

% february 2005 ---------------------------------------------
figure; prop = tfeb; u = ufeb; v = vfeb;
pcolor(lon(1,:),lat(:,1),prop);shading flat;hold on; cb = colorbar;
contour(xb,yb,zb,[-200 -1000],'w');
caxis(clim);
if quiv == 'y'; 
eval([QUIVER]); 
quiver(lonlim(1)+0.2,latlim(2)-0.2,1*sc,0*sc,0,'k');
text(lonlim(1)+0.2,latlim(2)-0.5,'1 m/s');
end 
plot(ncst(:,1),ncst(:,2),'k');
axis([lonlim latlim])
title(['hycom - temp - ',DEPTH,'m - feb/',YEAR,'']);
daspect([DASPECT])

% march 2005 ---------------------------------------------
figure; prop = tmar; u = umar; v = vmar;
pcolor(lon(1,:),lat(:,1),prop);shading flat;hold on; cb = colorbar;
contour(xb,yb,zb,[-200 -1000],'w');
caxis(clim);
if quiv == 'y'; 
eval([QUIVER]); 
quiver(lonlim(1)+0.2,latlim(2)-0.2,1*sc,0*sc,0,'k');
text(lonlim(1)+0.2,latlim(2)-0.5,'1 m/s');
end 
plot(ncst(:,1),ncst(:,2),'k');
axis([lonlim latlim])
title(['hycom - temp - ',DEPTH,'m - mar/',YEAR,'']);
daspect([DASPECT])

% february/march 2005 - OEII times  ---------------------------------------
figure; prop = toe2; u = uoe2; v = voe2;
pcolor(lon(1,:),lat(:,1),prop);shading flat;hold on; cb = colorbar;
contour(xb,yb,zb,[-200 -1000],'w');
caxis(clim);
if quiv == 'y'; 
eval([QUIVER]); 
quiver(lonlim(1)+0.2,latlim(2)-0.2,1*sc,0*sc,0,'k');
text(lonlim(1)+0.2,latlim(2)-0.5,'1 m/s');
end 
plot(ncst(:,1),ncst(:,2),'k');
axis([lonlim latlim])
title(['hycom - temp - ',DEPTH,'m - feb/mar/',YEAR,'']);
daspect([DASPECT])

if PRINT == 'y'; eval(['print -dpng ',figdir,'/',datatype,'/',id,'_',PROP,'_',DEPTH,'m_feb-mar_averages_',YEAR,'.png']);  
   tic; disp('Printing figure......');toc 
end




