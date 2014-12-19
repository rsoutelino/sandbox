%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT FOR VIZUALIZATION OF SST IMAGES  
% - HORIZONTAL MAPS
% - Rafael Soutelino - rsoutelino@gmail.com
% - last update: August, 2009
%
% - Most parameters to be changed are at 'SETTINGS' section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = 'proab1'
satelite = 'goes'; num = 30; 
%  satelite = 'modis'; num = 7;
pathname = ['/home/rafaelgs/doutorado/data/',satelite,'/',dataset]
lonlim = [-41 -34]; latlim = [-20 -9];
m_proj('mercator','lat',latlim,'lon',lonlim);
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');
load costa.mat;
load /usr/local/matlab/toolbox/lado/m_map/redblue;
% basic figures settings
PRINT = 'n'; % flag to print figures or not
DASPECT = [0.8 1 0.8];
DASPECT2 = [1 250 1];
s = 1; % scale size for vectors
ax = ['axis([min(xsec) max(xsec) -1000 0])'];

% END OF SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('==> .................. PLOTTING 8 DAY SST IMAGES ....................')
disp(' ')

sstm = [];

for k = 1:num

   if k < 10
   ncload([pathname,'/',satelite,'_sst8day_2007_000',num2str(k),'.nc'])
   else
   ncload([pathname,'/',satelite,'_sst8day_2007_00',num2str(k),'.nc'])
   end   

   disp(' ')
   disp('==> ............... MASKING CONTINENTAL SHELF .................')
   disp(' ')
   [lon,lat] = meshgrid(longitude-360,latitude);
   zb2 = griddata(xb,yb,zb,lon,lat);
   f = find(zb2 > -100); sst(f) = nan;
   sstm(k,:,:) =  sst ;

   h = figure('visible','off');
   pcolor(lon,lat,sst);shading flat;hold on;colorbar;caxis([23 29])
   plot(ncst(:,1),ncst(:,2),'k');
   contour(lon,lat,zb2,[-200 -1000],'k');
   axis([lonlim latlim])
   daspect([DASPECT])
   eval(['print -dpng ',pathname,'/figures/sst_',dataset,'_8day_',num2str(k),'.png'])

end


sstm = nanmean(sstm);
sstm = squeeze(sstm);

   figure
   pcolor(lon,lat,sstm);shading flat;hold on;colorbar;caxis([23 29])
   plot(ncst(:,1),ncst(:,2),'k');
   contour(lon,lat,zb2,[-200 -1000],'k');
   axis([lonlim latlim])
   daspect([DASPECT])
   eval(['print -dpng ',pathname,'/figures/sst_',dataset,'_mean.png'])



















