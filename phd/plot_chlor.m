%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT FOR VIZUALIZATION OF AVISO COMBINED MODIS-SEAWIFS 
%  COLOR IMAGES 
% - HORIZONTAL MAPS
% - Rafael Soutelino - rsoutelino@gmail.com
% - last update: August, 2009
%
% - Most parameters to be changed are at 'SETTINGS' section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = 'leste2'
pathname = ['/home/rafaelgs/doutorado/data/modis/',dataset]
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
disp('==> .................. PLOTTING 8 DAY COLOR IMAGES ....................')
disp(' ')

chlor = [];

for k = 1:7

   ncload([pathname,'/modis-seawifs_chlor8day_2005_000',num2str(k),'.nc'])
   
   disp(' ')
   disp('==> ............... MASKING CONTINENTAL SHELF .................')
   disp(' ')
   [lon,lat] = meshgrid(longitude,latitude);
   zb2 = griddata(xb,yb,zb,lon,lat);
   f = find(zb2 > -100); XSChlor(f) = nan;
   chlor(k,:,:) =  XSChlor;  

   h = figure('visible','off');
   pcolor(lon,lat,XSChlor);shading flat;hold on;colorbar;caxis([0 0.1])
   plot(ncst(:,1),ncst(:,2),'k');
   contour(lon,lat,zb2,[-200 -1000],'k');
   axis([lonlim latlim])
   daspect([DASPECT])
   eval(['print -dpng ',pathname,'/figures/chlor_',dataset,'_8day_',num2str(k),'.png'])

end


chlorm = nanmean(chlor);
chlorm = squeeze(chlorm);

   figure
   pcolor(lon,lat,chlorm);shading flat;hold on;colorbar;caxis([0 0.1])
   plot(ncst(:,1),ncst(:,2),'k');
   contour(lon,lat,zb2,[-200 -1000],'k');
   axis([lonlim latlim])
   daspect([DASPECT])
   eval(['print -dpng ',pathname,'/figures/chlor_',dataset,'_mean.png'])



















