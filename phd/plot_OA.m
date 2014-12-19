%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT FOR VIZUALIZATION OF HOPS OA NETCDF FILE OUTPUT
% - HORIZONTAL MAPS
% - VERTICAL SECTIONS
% - Rafael Soutelino - rsoutelino@gmail.com
% - last update: August-05, 2009
%
% - Most parameters to be changed are at 'SETTINGS' section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = 'leste2'
pathname = ['/home/rafaelgs/phd/data/',dataset,'/']
%  oafile = 'OA_leste2.nc'
%  oafile = 'OA_leste2abro2.nc'
%  oafile = 'OA_leste2abro2_levitus_roms.nc'
%  oafile = 'OA_levitus_marco_leste.nc'
oafile = 'OA_leste2_levitus.nc'
%  oafile = 'OA_leste2_levitus_cc50.nc'
%  oafile = 'OA_levitus_leste2.nc'
%  oafile = 'OA_proab1_levitus.nc'
%  oafile = 'OA_leste1_levitus.nc'
lonlim = [-45 -30]; latlim = [-24 -9];
z = [0:10:1000];
%  z = [0.,    10.,   20.,   30.,   50.,   75.,   100.,  125.,  150.,  200., ...
%       250.,  300.,  400.,  500.,  600.,  700.,  800.,  900.,  1000., 1100., ...
%       1200., 1300., 1400., 1500., 1750., 2000., 2500., 3000., 3500., 4000., ...
%       4500., 5000., 5500]; 
z = -z;
m_proj('mercator','lat',latlim,'lon',lonlim);
ref = 1000; % level of no motion stabilished at OA settings
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');
load costa.mat;
load /usr/local/matlab/toolbox/lado/m_map/redblue;
lT = [0 30];
lS = [33 38];
lV = [-1 1];

PLOT = menu('Choose the type of plot','Horizontal','Vertical');
if PLOT == 1;
   n = inputdlg('Pick a depth level [m]:   '); n = cell2mat(n); n = str2num(n);
   field = menu('Pick a field to plot: ','TEMP','SALT','PSI','TOPO');
   if field == 1; PROP = 'temp';
   elseif field == 2; PROP = 'salt';
   elseif field == 3; PROP = 'geopsi';
   else; PROP = 'topo';
   end
   INFO = ['xlabel(''Longitude''); ylabel(''Latitude''); title('' ',num2str(PROP),' - ',num2str(n),' m '')'];
else
   field = menu('Pick a field to plot: ','TEMP','SALT','VGEO');
   if field == 1; PROP = 'temp'; rg = lT;
   elseif field == 2; PROP = 'salt'; rg = lS;
   elseif field == 3; PROP = 'vgeo'; rg = lV; rg2 = -1:0.005:1;
   end
   sec = menu('What kind of section do you want to plot?','Zonal','Meridional');
   % ATTENTION!!  Must edit bellow lines if you are at northern hemisphere
   if sec == 1; 
      secc = inputdlg('Then pick the latitude of the section (S): '); 
      secc = cell2mat(secc); secc = str2num(secc); secc = -secc;
      INFO = ['xlabel(''Longitude''); ylabel(''Depth''); title('' ',num2str(PROP),' - ',num2str(-secc),' S '')'];
   else
      secc = inputdlg('Then pick the longitude of the section (W): ');
      secc = cell2mat(secc); secc = str2num(secc); secc = -secc;
      INFO = ['xlabel(''Latitude''); ylabel(''Depth''); title('' ',num2str(PROP),' - ',num2str(-secc),' W '')'];
   end
end

% basic figures settings
PRINT = 'y'; % flag to print figures or not
DASPECT = [1 1 1];
DASPECT2 = [1 250 1];
s = 1; % scale size for vectors
dx = 3; dy = 3; % subsampling fator for vectors
quiv = 's';
QUIVER = ['quiver(lon(1:dx:end,1:dy:end),lat(1:dx:end,1:dy:end),U(1:dx:end,1:dy:end)*s,V(1:dx:end,1:dy:end)*s,0,''k'');'];
secax = ['axis([min(xsec) max(xsec) -1000 0])'];
%  secax = ['axis([-39 -35 -1000 0])'];
%  secax = ['axis([-38 -34 -1000 0])'];
%  secax = ['axis([-22 -9 -1000 0])'];

% END OF SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('==> .................. LOADING FIELDS ....................')
disp(' ')
nc = netcdf([pathname,oafile]);
lon  = nc{'grid3'}(:,:,1);
lat  = nc{'grid3'}(:,:,2);
temp = nc{'temp'};  temp = temp(:);
salt = nc{'salt'};  salt = salt(:);
gpan = nc{'dynht'}; gpan = gpan(:); gpan = gpan*10;
[i,j,k] = size(temp);

for a = 1:k
   x(:,:,a) = lon(:,:);
   y(:,:,a) = lat(:,:);
end
[z2,x2] = meshgrid(z,lon(1,:));
for a = 1:i
   z3(a,:,:) = z2(:,:);
end

disp(' ')
disp('==> ...... COMPUTING STREAMFUNCTION AND GEOSTROPHIC VELOCITIES .......')
disp(' ')
disp(['    ==> Using ',num2str(ref),' dbar as reference .......'])

f0 = sw_f(mean(mean(lat)));
psi = gpan/f0;
ug=[];vg=[];
for a = 1:length(z);
   [ug(:,:,a),vg(:,:,a)] = psi2uv(lon,lat,psi(:,:,a));
end

disp(' ')
disp('==> .................. MASKING FIELDS ....................')
disp(' ')

h = griddata(xb,yb,zb,lon,lat);
for a = 1:length(z);
   f = find(h > z(a));
   temp2 = temp(:,:,a); temp2(f) = nan; temp(:,:,a) = temp2; clear temp2;
   salt2 = salt(:,:,a); salt2(f) = nan; salt(:,:,a) = salt2; clear salt2;
   psi2 = psi(:,:,a); psi2(f) = nan; psi(:,:,a) = psi2; clear psi2;
   ug2 = ug(:,:,a); ug2(f) = nan; ug(:,:,a) = ug2; clear ug2;
   vg2 = vg(:,:,a); vg2(f) = nan; vg(:,:,a) = vg2; clear vg2;
end

% =======================================================================

disp(' ')
disp('==> ............. 3D PLOT TO CHECK DATA SET ................')
disp(' ')
 
%  figure
%  p2 = patch(isocaps(x,y,z3,temp,-0.001), 'FaceColor', 'interp', 'EdgeColor', 'none');
%  p = patch(isosurface(x,y,z,temp,2));
%  hold on 
%  fill(ncst(:,1),ncst(:,2),[.7 .7 .7])  
%  contour(xb,yb,zb,[-200 -1000],'w')
%  daspect([.1 .1 30])
%  view(3)
%  colorbar('horiz')
%  grid on 
%  box off

if PLOT == 1 % ========================================================= 

disp(' ')
disp('==> .................. HORIZONTAL PLOTS ....................')
disp(' ')

n2 = near(z,-n,1);
T = squeeze(temp(:,:,n2)); S = squeeze(salt(:,:,n2)); PSI = squeeze(psi(:,:,n2));
U = squeeze(ug(:,:,n2)); V = squeeze(vg(:,:,n2));

if field == 1; prop = T; 
elseif field == 2; prop = S;
elseif field == 3; prop = PSI;
else; prop = h; 
end

figure
contourf(lon,lat,prop,30);shading flat;hold on; colorbar
contour(xb,yb,zb,[-200 -1000],'w');
if quiv == 's'; 
eval([QUIVER]); 
quiver(lonlim(1)+1,latlim(2)-1,1*s,0*s,0,'k');
text(lonlim(1)+1,latlim(2)-1.5,'1 m/s');
end 
plot(ncst(:,1),ncst(:,2),'k');
axis([lonlim latlim])
eval([INFO])
daspect([DASPECT])
if PRINT == 'y'; eval(['print -depsc figures/',PROP,'_',dataset,'_',num2str(n),'m.eps']);  end

else % ================================================================== 

disp(' ')
disp('==> .................. VERTICAL PLOTS ....................')
disp(' ')

if field == 1; prop = temp; 
elseif field == 2; prop = salt;
end

if sec == 1; 
   xsec = lat(:,1); fsec = near(xsec,secc,1);
   xsec = lon(1,:); xsec = reshape(xsec,1,length(xsec));
   if field == 3; prop = vg; end;
   prop = squeeze(prop(fsec,:,:)); prop = prop';
else; 
   xsec = lon(1,:); fsec = near(xsec,secc,1);
   xsec = lat(:,1); xsec = reshape(xsec,1,length(xsec));
   if field == 3; prop = ug; end; 
   prop = squeeze(prop(:,fsec,:)); prop = prop';
end

%======================================================================

figure
if field == 3;
   contourf(xsec,z,prop,rg2);shading flat;hold on;
   colormap(flipud(redblue));
else
   contourf(xsec,z,prop,50);shading flat;hold on;
end
caxis([rg]); colorbar
eval([secax]);
eval([INFO])
daspect([DASPECT2])
if PRINT == 'y'; 
   if sec == 1; 
      eval(['print -depsc figures/sec_',PROP,'_',dataset,'_',num2str(abs(secc)),'S.eps']);
   else   
      eval(['print -depsc figures/sec_',PROP,'_',dataset,'_',num2str(abs(secc)),'W.eps']);
   end
end

end


%========================================================================
% FURTHER ANALYSIS


disp(' ')
disp('==> ............. PLOTTING SIGMA LAYERS DEPTHS ................')
disp(' ')

%%% computing sigma_theta
%  dens = sw_dens0(salt,temp)-1000;
%  
%  iacas = [25.7 26.9]; % 
%  iaia = [26.9 27.53];
%  
%  acas = []; aia = [];
%  for ai = 1:i;
%     for aj = 1:j;
%        f1 = find(dens(ai,aj,:) > iacas(1) & dens(ai,aj,:) < iacas(end));
%        f2 = find(dens(ai,aj,:) > iaia(1) & dens(ai,aj,:) < iaia(end));
%        if isempty(f1) == 1
%           acas(ai,aj) = nan;
%        else
%           acas(ai,aj) = z(round(mean(f1)));
%        end
%        if isempty(f2) == 1
%           aia(ai,aj) = nan;
%        else
%           aia(ai,aj) = z(round(mean(f2)));
%        end
%     end
%  end
%  
%  figure
%  contourf(lon,lat,acas,50);shading flat;hold on; colorbar
%  contour(xb,yb,zb,[-200 -1000],'w');
%  if quiv == 's'; 
%  eval([QUIVER]); 
%  quiver(lonlim(1)+1,latlim(2)-1,1*s,0*s,0,'k');
%  text(lonlim(1)+1,latlim(2)-1.5,'1 m/s');
%  end 
%  plot(ncst(:,1),ncst(:,2),'k');
%  axis([lonlim latlim])
%  eval([INFO])
%  daspect([DASPECT])
%  
%  figure
%  contourf(lon,lat,aia,30);shading flat;hold on; colorbar
%  contour(xb,yb,zb,[-200 -1000],'w');
%  if quiv == 's'; 
%  eval([QUIVER]); 
%  quiver(lonlim(1)+1,latlim(2)-1,1*s,0*s,0,'k');
%  text(lonlim(1)+1,latlim(2)-1.5,'1 m/s');
%  end 
%  plot(ncst(:,1),ncst(:,2),'k');
%  axis([lonlim latlim])
%  eval([INFO])
%  daspect([DASPECT])

disp(' ')
disp('==> ....... PLOTTING MERIDIONAL SECTION OF MERIDIONAL VELOCITIES ................')
disp(' ')

vm = [];
for a = 1:k;
    vaux = vg(:,:,a);
    vaux = nanmean(vaux');
    vm = [vm; vaux];
end

figure
subplot(211)
contourf(xsec,z,vm,rg2);shading flat;hold on;
colormap(flipud(redblue));
caxis([rg/10]); colorbar
eval([secax]);
eval([INFO])
daspect([DASPECT2])

disp(' ')
disp('==> ....... MERIDIONAL TRANSPORT COMPUTATION ................')
disp(' ')

zmin = 0;    % maximum level to compute transport
zmax = 1000; % maximum level to compute transport
nmin = near(z,-zmin,1);
nmax = near(z,-zmax,1);

dx = mean(sw_dist(lat(1,:),lon(1,:),'km'))*1000;
dz = abs(mean(diff(z)));
tm = [];
for a = 1:i
   vaux = squeeze(vg(a,:,nmin:nmax));
   vaux = nansum(nansum((vaux .* dx*dz))) * 1e-6; vaux = round(vaux*10)/10;
   tm = [tm vaux];
end

subplot(212)
plot(xsec,tm,'k'); hold on; grid on;
plot([xsec(1) xsec(end)],[0 0],'k--');
set(gca,'xlim',[xsec(1) xsec(end)])

if PRINT == 'y'; 
   if sec == 1; 
      eval(['print -depsc figures/secbif_',PROP,'_',dataset,'_',num2str(abs(secc)),'S.eps']);
   else   
      eval(['print -depsc figures/secbif_',PROP,'_',dataset,'_',num2str(abs(secc)),'W.eps']);
   end
end


%  end














