%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting HYCOM GLOBAL outputs
% Rafael Soutelino - May 2010
% rsoutelino@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYCOM + NCODA Global 1/12 Analysis (expt_60.5) (assimilates data) Nov-2003 to Dec-2006
url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5'; 

ltemp = 3:0.04:29;   % colorbar scale for temperature
lsalt = 30:0.01:37;   % colorbar scale for salinity
leta=-0.2:0.01:0.2;   % colorbar scale for surface elevation
lu = -1:0.01:1;       % colorbar scale for zonal velocity
lv = -1:0.01:1; % colorbar scale for meridional velocity
lw = -0.001:0.00001:0.001; % colorbar scale for vertical velocity
lT = [0 30];
lS = [33 38];
lV = [-0.5 0.5];  
% settings to interpolate s-coord to z-coord  
dz = 10;
isob=[-200 -1000];    % isobaths to be ploted
lonlim = [-41 -34]; latlim = [-22 -12];
%  lonlim = [-41 -34]; latlim = [-21 -10];
m_proj('mercator','long',lonlim,'lat',latlim,'on');
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');
load costa.mat           % m_map inicialization
load /Users/rsoutelino/rsoutelino/mymatlab/lado/m_map/redblue
load /Users/rsoutelino/rsoutelino/misc/levitus_levels.txt;
woa = levitus_levels(:,2); clear levitus_levels;
figdir = 'figures'   % directory name to save the figures
datadir = 'outputs'
depth = 50; 
id = 'hycom'          % run ID
datatype = 'hycom_glb08_exp60.5'      % output type: HIS or AVG?
YEAR = '2005';
MONTH = 'mar'
day = [1:30]; % input('Pick a day:   ');  


% basic figures settings
PRINT = 'y'; % flag to print figures or not
DASPECT = [1 1 1];
DASPECT2 = [1 250 1];
sc = 1; % scale size for vectors
quiv = 'y';
dx = 2; dy = 2;      % vectors ploting interval
QUIVER = ['quiver(lon(1:dx:end,1:dy:end),lat(1:dx:end,1:dy:end),u(1:dx:end,1:dy:end)*sc,v(1:dx:end,1:dy:end)*sc,0,''k'');'];
secax = ['axis([min(xsec) max(xsec) -1000 0])'];

% velocity over prop plots
PROP = 'temp'; 
prop2 = ['prop = T;'];
clim = [15 28]; 
INFO = ['xlabel(''Longitude''); ylabel(''Latitude''); title(['' ',id,' - ',PROP,' - '',DEPTH,''m - ',MONTH,'/'',DAY,''/',YEAR,' ''])'];

proc = 'y'; % flag to perform or not further analysis

% saving variables for further mean calculations
if proc == 'y' ucum = []; vcum = []; tcum = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading and preparing the output variables --------------------------------------

for k = day
DAY = num2str(k);

disp('   ')
disp('Reading the output netcdf file.......')
disp('   ')

tic
ofile  = [datadir,'/',id,'_',[DAY MONTH YEAR],'.nc'];
nco  = netcdf(ofile);
T    =  nco{'TEMPERATURE'}(:,:,:);
S    =  nco{'SALINITY'}(:,:,:);
U    =  nco{'U'}(:,:,:);
V    =  nco{'V'}(:,:,:);
lon    =  nco{'LONGITUDE'}(:,:,:); lon = lon-360;
lat    =  nco{'LATITUDE'}(:,:,:);	
%  ubt  =  nco{'ubar'}(day,:,:);
%  vbt  =  nco{'vbar'}(day,:,:);
%  eta  =  nco{'zeta'}(day,:,:);
toc

f = find(T > 50);
T(f) = nan; S(f) = nan; U(f) = nan; V(f) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing horizontal fields
f = near(woa,depth,1); DEPTH = num2str(woa(f));
eval([prop2]);
prop = squeeze(prop(f,:,:)); u = squeeze(U(f,:,:)); v = squeeze(V(f,:,:));  

%%% plotting horizontal fields -=========================================================
disp(' ')
disp('==> .................. HORIZONTAL PLOTS ....................')
disp(' ')

if length(day) > 1; 
   h = figure('visible','off'); 
else
   figure
end

pcolor(lon(1,:),lat(:,1),prop);shading flat;hold on; colorbar
contour(xb,yb,zb,[-200 -1000],'w');
caxis(clim);
if quiv == 'y'; 
eval([QUIVER]); 
quiver(lonlim(1)+0.2,latlim(2)-0.2,1*sc,0*sc,0,'k');
text(lonlim(1)+0.2,latlim(2)-0.5,'1 m/s');
end 
plot(ncst(:,1),ncst(:,2),'k');
axis([lonlim latlim])
eval([INFO])
daspect([DASPECT])
if PRINT == 'y'; eval(['print -dpng ',figdir,'/',datatype,'/',id,'_',PROP,'_',DEPTH,'m_',MONTH,DAY,'_',YEAR,'.png']);  
   tic; disp('Printing figure......');toc 
end
if length(day) > 1; close all ; end

% saving variables for further mean calculations
if proc == 'y'; tcum(k,:,:) = prop; ucum(k,:,:) = u; vcum(k,:,:) = v; end

end

% saving variables for further mean calculations
if proc == 'y'; eval(['save ',datatype,'_',MONTH,YEAR,'_',DEPTH,'m.mat']); 
   tic; disp('Saving mat file.........'); toc;
end



























