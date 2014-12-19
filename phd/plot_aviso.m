%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT FOR VIZUALIZATION OF AVISO ALTIMETRY 
%      PRODUCT AND COMPUTING MEAN FIELDS
% - HORIZONTAL MAPS
% - Rafael Soutelino - rsoutelino@gmail.com
% - last update: August, 2009
%
% - Most parameters to be changed are at 'SETTINGS' section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = 'proab1'

if dataset == 'leste2' | dataset == 'brazil';
   matfile = 'absvel_2005.mat'
   mm = [2 3];
elseif dataset == 'proab1'
   matfile = 'absvel_2007.mat' 
   mm = [9];
end

pathname = ['/home/rafaelgs/doutorado/data/aviso/',dataset,'/matfiles/']
load([pathname,matfile]);
lonlim = [-45 -30]; latlim = [-24 -9];
m_proj('mercator','lat',latlim,'lon',lonlim);
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');
load costa.mat;
load /usr/local/matlab/toolbox/lado/m_map/redblue;
% basic figures settings
PRINT = 'n'; % flag to print figures or not
DASPECT = [1 1 1];
DASPECT2 = [1 250 1];
s = 1; % scale size for vectors
dx = 3; dy = 3; % subsampling fator for vectors
ax = ['axis([min(xsec) max(xsec) -1000 0])'];

% END OF SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('==> .................. COMPUTING MEAN FIELD ....................')
disp(' ')

day = datestr(day,6);
mon = day(:,1:2); mon = str2num(mon);

f = [];
for k = 1:length(mm)
   F = find(mon == mm(k));
   f = [F;f];
end

um = squeeze(nanmean(u(f,:,:))); um = um';
vm = squeeze(nanmean(v(f,:,:))); vm = vm';


lon2 = lon(1):0.1:lon(end);
lat2 = lat(1):0.1:lat(end);
[lon,lat] = meshgrid(lon,lat);
[lon2,lat2] = meshgrid(lon2,lat2);

um = griddata(lon,lat,um,lon2,lat2);
vm = griddata(lon,lat,vm,lon2,lat2);
mag = sqrt(um.^2 + vm.^2);

disp(' ')
disp('==> ............... MASKING CONTINENTAL SHELF .................')
disp(' ')

zb = griddata(xb,yb,zb,lon2,lat2);
f = find(zb > -1000);
um(f) = nan; vm(f) = nan; mag(f) = nan;

pcolor(lon2,lat2,mag); shading flat; colorbar; hold on
colormap(redblue(32:end,:));
str = ['(1:dx:end,1:dy:end)'];
eval(['quiver(lon2' str ',lat2' str ',um' str ',vm' str ',3,''k'')']);
plot(ncst(:,1),ncst(:,2),'k');
contour(lon2,lat2,zb,[-200 -1000],'k');
caxis([0 50])
axis([lonlim latlim])
daspect([DASPECT])




















