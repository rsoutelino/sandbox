function mat2roms_rps(theMatFile, theROMSFile)
% mat2roms_rps -- Convert from Mat-file to ROMS file.
%  mat2roms('theMatFile', 'theROMSFile') converts
%   'theMatFile' (mostly RHO information) to 'theROMSFile'.
%
% NOTE: This is a modified version of Chuck Denham's mat2roms function.
% In this version the only requirements for the input mat file are: 
% s.rho.lon
% s.rho.lat
% s.rho.depth,  where (depth == nan) indicates land points
%
% The projected coordinate is calculated via the projection
% specified in this routine. This is *only* used to interpolate
% the grid to find the lon/lat values of the u,v and psi
% points and to write out x_rho, y_rho, x_u, y_u, etc, which
% are not used by ROMS, but could be useful for
% plotting.   The only important quantities for ROMS are the
% metric factors pm and pn and the angles, which 
% are calculated from the lon/lat values
% assuming a spherical earth, not from distances in the 
% projected coordinates.
%
% -Rich Signell (21-May-2003) rsignell@usgs.gov
%
% This routine calls m_proj from the m_map toolkit and
% also sw_dist and sw_f from the seawater toolkit (both
% available at http://sea-mat.whoi.edu
 
% Copyright (C) 2002 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 22-May-2002 16:25:10.
% Updated    17-Jun-2002 16:49:32.
% Updated by RPS 21-May-2003

isMacintosh = ~isunix & any(findstr(lower(computer), 'mac'));

if nargin < 1, theMatFile = '*.mat'; end
if nargin < 2, theROMSFile = 'roms_model_grid.nc'; end

% Get the file names.

if any(theMatFile == '*')
	help(mfilename)
	theFilterSpec = theMatFile;
	thePrompt = 'Select a Mat-File';
	[theFile, thePath] = uigetfile(theFilterSpec, thePrompt);
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theMatFile = [thePath theFile];
end

if any(theROMSFile == '*')
	theFilterSpec = theROMSFile;
	thePrompt = 'Save As ROMS File';
	[theFile, thePath] = uiputfile(theFilterSpec, thePrompt);
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theROMSFile = [thePath theFile];
end

if isequal(theMatFile, theROMSFile)
	disp([' ## Must not select same file for input and output.'])
	return
end

% Load the Mat-File.

s = load(theMatFile);
if isempty(s)
	disp([' ## Mat-File is empty.'])
	return
end

% We need to pay attention to the (i, j) coordinates
%  stored in the file.  They are our only clue on the
%  actual orientation of the data.  The "i" index
%  goes left-to-right in the grid-world, while
%  the "j" index runs bottom-to-top.

% WetCDF on.

if isMacintosh
	eval('wetcdf on')
end

% Open the ROMS File.

nc = netcdf(theROMSFile, 'clobber');
if isempty(nc)
	disp([' ## Unable to open ROMS NetCDF output file.'])
	return
end

% Populate the ROMS File.
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
 
nc.type = ncchar('Gridpak file');
theGridTitle = 'ROMS Model Grid';
nc.gridid = theGridTitle;
nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);

nc.CPP_options = ncchar('DCOMPLEX, DBLEPREC, NCARG_32, PLOTS,');
name(nc.CPP_options, 'CPP-options')

% Dimensions:


% mapping

[m, n] = size(s.rho.lon);

% The xi direction (left-right):

LP = n;   % The rho dimension.
L = LP-1; % The psi dimension.

% The eta direction (up-down):

MP = m;   % The rho dimension.
M = MP-1; % The psi dimension.

disp(' ## Defining Dimensions...')
 
nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;

nc('one') = 1;
nc('two') = 2;
nc('bath') = 0; %% (record dimension)
 
%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')
 
nc{'xl'} = ncdouble('one'); %% 1 element.
nc{'xl'}.long_name = ncchar('domain length in the XI-direction');
nc{'xl'}.units = ncchar('meter');
 
nc{'el'} = ncdouble('one'); %% 1 element.
nc{'el'}.long_name = ncchar('domain length in the ETA-direction');
nc{'el'}.units = ncchar('meter');
 
nc{'JPRJ'} = ncchar('two'); %% 2 elements.
nc{'JPRJ'}.long_name = ncchar('Map projection type');

nc{'JPRJ'}.option_ME_ = ncchar('Mercator');
nc{'JPRJ'}.option_ST_ = ncchar('Stereographic');
nc{'JPRJ'}.option_LC_ = ncchar('Lambert conformal conic');
name(nc{'JPRJ'}.option_ME_, 'option(ME)')
name(nc{'JPRJ'}.option_ST_, 'option(ST)')
name(nc{'JPRJ'}.option_LC_, 'option(LC)')
 
nc{'PLAT'} = ncfloat('two'); %% 2 elements.
nc{'PLAT'}.long_name = ncchar('Reference latitude(s) for map projection');
nc{'PLAT'}.units = ncchar('degree_north');
 
nc{'PLONG'} = ncfloat('one'); %% 1 element.
nc{'PLONG'}.long_name = ncchar('Reference longitude for map projection');
nc{'PLONG'}.units = ncchar('degree_east');
 
nc{'ROTA'} = ncfloat('one'); %% 1 element.
nc{'ROTA'}.long_name = ncchar('Rotation angle for map projection');
nc{'ROTA'}.units = ncchar('degree');
 
nc{'JLTS'} = ncchar('two'); %% 2 elements.
nc{'JLTS'}.long_name = ncchar('How limits of map are chosen');

nc{'JLTS'}.option_CO_ = ncchar('P1, .. P4 define two opposite corners ');
nc{'JLTS'}.option_MA_ = ncchar('Maximum (whole world)');
nc{'JLTS'}.option_AN_ = ncchar('Angles - P1..P4 define angles to edge of domain');
nc{'JLTS'}.option_LI_ = ncchar('Limits - P1..P4 define limits in u,v space');
name(nc{'JLTS'}.option_CO_, 'option(CO)')
name(nc{'JLTS'}.option_MA_, 'option(MA)')
name(nc{'JLTS'}.option_AN_, 'option(AN)')
name(nc{'JLTS'}.option_LI_, 'option(LI)')
 
nc{'P1'} = ncfloat('one'); %% 1 element.
nc{'P1'}.long_name = ncchar('Map limit parameter number 1');
 
nc{'P2'} = ncfloat('one'); %% 1 element.
nc{'P2'}.long_name = ncchar('Map limit parameter number 2');
 
nc{'P3'} = ncfloat('one'); %% 1 element.
nc{'P3'}.long_name = ncchar('Map limit parameter number 3');
 
nc{'P4'} = ncfloat('one'); %% 1 element.
nc{'P4'}.long_name = ncchar('Map limit parameter number 4');
 
nc{'XOFF'} = ncfloat('one'); %% 1 element.
nc{'XOFF'}.long_name = ncchar('Offset in x direction');
nc{'XOFF'}.units = ncchar('meter');
 
nc{'YOFF'} = ncfloat('one'); %% 1 element.
nc{'YOFF'}.long_name = ncchar('Offset in y direction');
nc{'YOFF'}.units = ncchar('meter');
 
nc{'depthmin'} = ncshort('one'); %% 1 element.
nc{'depthmin'}.long_name = ncchar('Shallow bathymetry clipping depth');
nc{'depthmin'}.units = ncchar('meter');
 
nc{'depthmax'} = ncshort('one'); %% 1 element.
nc{'depthmax'}.long_name = ncchar('Deep bathymetry clipping depth');
nc{'depthmax'}.units = ncchar('meter');
 
nc{'spherical'} = ncchar('one'); %% 1 element.
nc{'spherical'}.long_name = ncchar('Grid type logical switch');
nc{'spherical'}.option_T_ = ncchar('spherical');
nc{'spherical'}.option_F_ = ncchar('Cartesian');
name(nc{'spherical'}.option_T_, 'option(T)')
name(nc{'spherical'}.option_F_, 'option(F)')
 
nc{'hraw'} = ncdouble('bath', 'eta_rho', 'xi_rho'); %% 0 elements.
nc{'hraw'}.long_name = ncchar('Working bathymetry at RHO-points');
nc{'hraw'}.units = ncchar('meter');
nc{'hraw'}.field = ncchar('bath, scalar');
 
nc{'h'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'h'}.long_name = ncchar('Final bathymetry at RHO-points');
nc{'h'}.units = ncchar('meter');
nc{'h'}.field = ncchar('bath, scalar');
 
nc{'f'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'f'}.long_name = ncchar('Coriolis parameter at RHO-points');
nc{'f'}.units = ncchar('second-1');
nc{'f'}.field = ncchar('Coriolis, scalar');
 
nc{'pm'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
nc{'pm'}.units = ncchar('meter-1');
nc{'pm'}.field = ncchar('pm, scalar');
 
nc{'pn'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
nc{'pn'}.units = ncchar('meter-1');
nc{'pn'}.field = ncchar('pn, scalar');
 
nc{'dndx'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dndx'}.long_name = ncchar('xi derivative of inverse metric factor pn');
nc{'dndx'}.units = ncchar('meter');
nc{'dndx'}.field = ncchar('dndx, scalar');
 
nc{'dmde'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dmde'}.long_name = ncchar('eta derivative of inverse metric factor pm');
nc{'dmde'}.units = ncchar('meter');
nc{'dmde'}.field = ncchar('dmde, scalar');
 
nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'x_rho'}.long_name = ncchar('x location of RHO-points');
nc{'x_rho'}.units = ncchar('meter');
 
nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'y_rho'}.long_name = ncchar('y location of RHO-points');
nc{'y_rho'}.units = ncchar('meter');
 
nc{'x_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'x_psi'}.long_name = ncchar('x location of PSI-points');
nc{'x_psi'}.units = ncchar('meter');
 
nc{'y_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'y_psi'}.long_name = ncchar('y location of PSI-points');
nc{'y_psi'}.units = ncchar('meter');
 
nc{'x_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'x_u'}.long_name = ncchar('x location of U-points');
nc{'x_u'}.units = ncchar('meter');
 
nc{'y_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'y_u'}.long_name = ncchar('y location of U-points');
nc{'y_u'}.units = ncchar('meter');
 
nc{'x_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'x_v'}.long_name = ncchar('x location of V-points');
nc{'x_v'}.units = ncchar('meter');
 
nc{'y_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'y_v'}.long_name = ncchar('y location of V-points');
nc{'y_v'}.units = ncchar('meter');
 
nc{'lat_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lat_rho'}.long_name = ncchar('latitude of RHO-points');
nc{'lat_rho'}.units = ncchar('degree_north');
 
nc{'lon_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lon_rho'}.long_name = ncchar('longitude of RHO-points');
nc{'lon_rho'}.units = ncchar('degree_east');
 
nc{'lat_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lat_psi'}.long_name = ncchar('latitude of PSI-points');
nc{'lat_psi'}.units = ncchar('degree_north');
 
nc{'lon_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lon_psi'}.long_name = ncchar('longitude of PSI-points');
nc{'lon_psi'}.units = ncchar('degree_east');
 
nc{'lat_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lat_u'}.long_name = ncchar('latitude of U-points');
nc{'lat_u'}.units = ncchar('degree_north');
 
nc{'lon_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lon_u'}.long_name = ncchar('longitude of U-points');
nc{'lon_u'}.units = ncchar('degree_east');
 
nc{'lat_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lat_v'}.long_name = ncchar('latitude of V-points');
nc{'lat_v'}.units = ncchar('degree_north');
 
nc{'lon_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lon_v'}.long_name = ncchar('longitude of V-points');
nc{'lon_v'}.units = ncchar('degree_east');
 
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nc{'mask_rho'}.option_0_ = ncchar('land');
nc{'mask_rho'}.option_1_ = ncchar('water');
name(nc{'mask_rho'}.option_0_, 'option(0)')
name(nc{'mask_rho'}.option_1_, 'option(1)')
 
nc{'mask_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'mask_u'}.long_name = ncchar('mask on U-points');
nc{'mask_u'}.option_0_ = ncchar('land');
nc{'mask_u'}.option_1_ = ncchar('water');
name(nc{'mask_u'}.option_0_, 'option(0)')
name(nc{'mask_u'}.option_1_, 'option(1)')
%		 nc{'mask_u'}.FillValue_ = ncdouble(1);
 
nc{'mask_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'mask_v'}.long_name = ncchar('mask on V-points');
nc{'mask_v'}.option_0_ = ncchar('land');
nc{'mask_v'}.option_1_ = ncchar('water');
name(nc{'mask_v'}.option_0_, 'option(0)')
name(nc{'mask_v'}.option_1_, 'option(1)')
%		 nc{'mask_v'}.FillValue_ = ncdouble(1);
 
nc{'mask_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'mask_psi'}.long_name = ncchar('mask on PSI-points');
nc{'mask_psi'}.option_0_ = ncchar('land');
nc{'mask_psi'}.option_1_ = ncchar('water');
name(nc{'mask_psi'}.option_0_, 'option(0)')
name(nc{'mask_psi'}.option_1_, 'option(1)')
%		 nc{'mask_psi'}.FillValue_ = ncdouble(1);

% Now, what about depths: "h" and "hraw".  <== DEPTHS.

nc{'angle'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'angle'}.long_name = ncchar('angle between xi axis and east');
nc{'angle'}.units = ncchar('radians');

% Variables.

% Fill the variables with data.

disp(' ## Filling Variables...')

if ~isfield(s, 'projection') | isempty(s.projection)
	s.projection = 'Mercator';
end
projection = s.projection;

switch lower(projection)
case 'mercator'
	theProjection = 'ME';
case 'stereographic'
	theProjection = 'ST';
case 'lambert conformal conic'
	theProjection = 'LC';
otherwise
	theProjection = '??';
end

nc{'JPRJ'}(1:2) = theProjection;
nc{'spherical'}(:) = 'T';   % T or F -- uppercase okay?

% Set projection parameters
%stdlt1o =   60.00
%stdlt2o =   30.00
%stdlono =   80.00
%m_proj('Lambert Conformal Conic','clon',stdlono,'para',[stdlt1o stdlt2o]);

m_proj(projection);

% Cartesian coordinates.
lon=s.rho.lon;
lat=s.rho.lat;

[x,y] = m_ll2xy(lon,lat);
radius=6371229.0;  % earth radius used by COAMPS
x=x*radius;
y=y*radius;

grid_x = interp2(x, 1);
grid_y = interp2(y, 1);

% interpolate Cartesian coordinates to get Geographic Coordinates at
% non-grid cell centers

[geogrid_lon, geogrid_lat] = m_xy2ll(grid_x/radius, grid_y/radius);   % Degrees.

xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));
nc{'xl'}(:) = xl;
nc{'el'}(:) = el;

f = sw_f(lat);
nc{'f'}(:) = f;

% Handy indices.

[m2, n2] = size(grid_x);

i_rho = 1:2:n2;
j_rho = 1:2:m2;
i_psi = 2:2:n2;
j_psi = 2:2:m2;
i_u = 2:2:n2;
j_u = 1:2:m2;
i_v = 1:2:n2;
j_v = 2:2:m2;

% Locations.

nc{'x_rho'}(:) = x;
nc{'y_rho'}(:) = y;

nc{'x_psi'}(:) = grid_x(j_psi, i_psi);
nc{'y_psi'}(:) = grid_y(j_psi, i_psi);

nc{'x_u'}(:) = grid_x(j_u, i_u);
nc{'y_u'}(:) = grid_y(j_u, i_u);

nc{'x_v'}(:) = grid_x(j_v, i_v);
nc{'y_v'}(:) = grid_y(j_v, i_v);

nc{'lon_rho'}(:) = geogrid_lon(j_rho, i_rho);
nc{'lat_rho'}(:) = geogrid_lat(j_rho, i_rho);

nc{'lon_psi'}(:) = geogrid_lon(j_psi, i_psi);
nc{'lat_psi'}(:) = geogrid_lat(j_psi, i_psi);

lon_u = geogrid_lon(j_u, i_u);
lat_u = geogrid_lat(j_u, i_u);
nc{'lon_u'}(:) = lon_u;
nc{'lat_u'}(:) = lat_u;

lon_v = geogrid_lon(j_v, i_v);
lat_v = geogrid_lat(j_v, i_v);
nc{'lon_v'}(:) = lon_v;
nc{'lat_v'}(:) = lat_v;

% Metric factors: use ones if absent.

if ~isfield(s.rho, 'dx') | isempty(s.rho.dx)
	s.rho.dx = ones(size(nc{'pm'}));
end
if ~isfield(s.rho, 'dy') | isempty(s.rho.dy)
	s.rho.dy = ones(size(nc{'pn'}));
end

% dx = s.rho.dx;    % NO!  This is wrong.  Should calculate
% dy = s.rho.dy;    % dx and dy on spherical earth, not in projected coords

% use sw_dist to compute metric distances between lon/lon points

lat_u=lat_u.';
lon_u=lon_u.';
[dx,ang]=sw_dist(lat_u(:),lon_u(:),'km');

dx=[dx(:); dx(end)]*1000;  % km==> m
dx=reshape(dx,n-1,m);
dx=dx.';
dx=[dx(:,1) dx(:,(1:(end-1))) dx(:,end-1)];

ang=[ang(:); ang(end)];
ang=reshape(ang,n-1,m);
ang=ang.';
ang=[ang(:,1) ang(:,(1:(end-1))) ang(:,end-1)];

dy=sw_dist(lat_v(:),lon_v(:),'km');
dy=[dy(:); dy(end)]*1000;   % km ==> m 
dy=reshape(dy,m-1,n);
dy=[dy(1,:); dy((1:(end-1)),:); dy(end-1,:)];

dx(dx == 0) = NaN;   % Shouldn't be any zeros here.
dy(dy == 0) = NaN;

pm = 1 ./ dx;   % Note reciprocals.
pn = 1 ./ dy;

% Does ROMS tolerate NaNs and/or Infs?
%  Do we have to substitute a NetCDF "_FillValue"?

nc{'pm'}(:) = pm;
nc{'pn'}(:) = pn;

dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));

nc{'dmde'}(:) = dmde;
nc{'dndx'}(:) = dndx;


% Depths at RHO points.

bathymetry = s.rho.depth;
mask = isnan(bathymetry) | (bathymetry == -99999);   % ECOM land = -99999.

bathymetry(mask) = min(bathymetry(:));   % set bathy under land to min depth
water = 1 - mask;   % ROMS requires water = 1, land= 0.
if ~isempty(bathymetry)
	nc{'h'}(:) = bathymetry;   % ROMS depths are positive.
end

% Angles at RHO points
nc{'angle'}(:) = ang*pi/180;   % Radians.

% The u, v, and psi masks rely directly on
%  the mask of the surrounding rho points.

nc{'mask_rho'}(:) = water;
nc{'mask_u'}(:) = water(:, 1:end-1) & water(:, 2:end);
nc{'mask_v'}(:) = water(1:end-1, :) & water(2:end, :);
nc{'mask_psi'}(:) = water(1:end-1, 1:end-1) &...
		water(1:end-1, 2:end) & ...
     		water(2:end, 1:end-1) & ...
  		water(2:end, 2:end);

% Close the ROMS File.

if ~isempty(close(nc))
	disp(' ## Unable too close the ROMS output file.')
end

% WetCDF off.

if isMacintosh
	eval('wetcdf off')
end

% ---------- romsnan ---------- %

function r = romsnan

% romsnan -- NaN value prefered by ROMS.
%  romsnan (no argument) returns the value prefered
%   by ROMS in place of NaN.  *** I believe it is
%   the NetCDF _FillValue for doubles.  (CRD)

f = ncfillvalues;
r = f.ncdouble;
