#!/usr/bin/env python
#
# Creates initial conditions netCDF file for ROMS
#
# Rafael Soutelino - rsoutelino@gmail.com
#
# Using some material from Matlab scripts by 
# "Copyright (c) 2003 UCLA - Patrick Marchesiello"
#
# Last modification: Aug, 2010
#####################################################################

print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 
# IMPORTING MODULES #################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import delaunay
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
import scipy.io as sp

# classes and functions to the computings
from roms_setup import run_setup, zlev, ztosigma 

#####################################################################

# SCRIPT START ######################################################

# Basic Settings:

filenamestr = '_ini.nc'
filetypestr = 'ROMS Initial Conditions file'

# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# OA-created netcdf initial T, S file 
# grid netcdf file

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
run = run_setup('../phd5_run.setup')

print ' \n' + '==> ' + '  READING OA NETCDF FILE ...\n' + ' '
datafile = sp.loadmat(run.datadir + run.ini_filename)

# assigning some variables from data file
Zlev   = datafile['depth'][:].ravel() * (-1)
Zlev[0] = 0
N1     = Zlev.size
lont   = datafile['lont'][:]
latt   = datafile['latt'][:]
lonu   = datafile['lonu'][:]
latu   = datafile['latu'][:]
u      = datafile['U'][:]
v      = datafile['V'][:]
zeta   = datafile['zeta'][:]
temp   = datafile['temp'][:]
salt   = datafile['salt'][:]

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(run.rundir + run.run_name + '_grd.nc')

# assigning some variables from grid file
rlon2  = grdfile.variables['lon_rho'][:]
rlat2  = grdfile.variables['lat_rho'][:]
vlon2  = grdfile.variables['lon_v'][:]
vlat2  = grdfile.variables['lat_v'][:]
ulon2  = grdfile.variables['lon_u'][:]
ulat2  = grdfile.variables['lat_u'][:]
angle  = grdfile.variables['angle'][:]
h2     = grdfile.variables['h'][:]
rmask2 = grdfile.variables['mask_rho'][:]

# DOING COMPUTATIONS TO INTERPOLATE THE FIELDS TO ROMS GRID #########

# Modify the bathymetry
f     = np.where(h2 >= 5000)
h2[f] = 5000; del f

Jrho, Irho = rlon2.shape
Mr2, Lr2   = rlon2.shape
Lu2 = Lr2-1; Mu2 = Mr2
Lv2 = Lr2;   Mv2 = Mr2-1
cosa = np.cos(angle); sina = np.sin(angle); del angle
rmask2 = np.ma.masked_where(rmask2 == 0, rmask2)

hu = griddata(rlon2.ravel(), rlat2.ravel(), h2.ravel(), ulon2, ulat2)
hv = griddata(rlon2.ravel(), rlat2.ravel(), h2.ravel(), vlon2, vlat2)
[Zsig,dZsig]   = zlev(h2,run.theta_s,run.theta_b,run.tcline,run.klevels)
[ZsigU,dZsigU] = zlev(hu,run.theta_s,run.theta_b,run.tcline,run.klevels)
[ZsigV,dZsigV] = zlev(hv,run.theta_s,run.theta_b,run.tcline,run.klevels)

### Interpolating T, S to ROMS 3D S-COORD grid ###############################
# Currently doing only TEMP and SALT. ZETA, UBAR, VBAR, U, V are zero arrays
t = np.arange(0, run.time);

lN   = run.klevels
lt   = np.size(run.time)
UBAR = np.zeros([lt, Mu2, Lu2])
VBAR = np.zeros([lt, Mv2, Lv2])

print ' \n' + '==> ' + '  INTERPOLATING FREE-SURFACE ...\n' + ' '
ZETA = griddata(lont.ravel(),latt.ravel(),zeta.ravel(),rlon2,rlat2)
ZETA.shape = (1, Jrho, Irho)

z2    = np.zeros([N1, Jrho, Irho])
Zlev2 = np.zeros([N1, 1])

print ' \n' + '==> ' + '  INTERPOLATING TEMPERATURE ...\n' + ' '
for k in np.arange(0, N1, 1): 
    print 'TEMP:  Z Level = ' + str(-1*Zlev[k]) + ' m'
    z1  = np.squeeze(temp[k,:,:])
    z2[N1-k-1,:,:] = griddata(lont.ravel(),latt.ravel(),z1.ravel(),rlon2,rlat2)
    Zlev2[N1-k-1]  = Zlev[k]

print ' \n' + '==> ' + '  INTERPOLATING TEMP FROM Z --> S COORD ...\n' + ' '
TEMP = ztosigma(z2,Zsig,Zlev2); del z1, z2

###

z2    = np.zeros([N1, Jrho, Irho])
print ' \n' + '==> ' + '  INTERPOLATING SALINITY ...\n' + ' '
for k in np.arange(0, N1, 1): 
    print 'SALT:  Z Level = ' + str(-1*Zlev[k]) + ' m'
    z1  = np.squeeze(salt[k,:,:])
    z2[N1-k-1,:,:] = griddata(lont.ravel(),latt.ravel(),z1.ravel(),rlon2,rlat2)
    Zlev2[N1-k-1]  = Zlev[k]

print ' \n' + '==> ' + '  INTERPOLATING SALT FROM Z --> S COORD ...\n' + ' '
SALT = ztosigma(z2,Zsig,Zlev2);

###

z2    = np.zeros([N1, Mu2, Lu2])
print ' \n' + '==> ' + '  INTERPOLATING U-velocity ...\n' + ' '
for k in np.arange(0, N1, 1): 
    print 'U-Vel:  Z Level = ' + str(-1*Zlev[k]) + ' m'
    z1  = np.squeeze(u[k,:,:])
    z2[N1-k-1,:,:] = griddata(lonu.ravel(),latu.ravel(),z1.ravel(),ulon2,ulat2)
    Zlev2[N1-k-1]  = Zlev[k]

print ' \n' + '==> ' + '  INTERPOLATING V-Vel FROM Z --> S COORD ...\n' + ' '
U   = ztosigma(z2,ZsigU,Zlev2);

###

z2    = np.zeros([N1, Mv2, Lv2])
print ' \n' + '==> ' + '  INTERPOLATING V-velocity ...\n' + ' '
for k in np.arange(0, N1, 1): 
    print 'V-Vel:  Z Level = ' + str(-1*Zlev[k]) + ' m'
    z1  = np.squeeze(v[k,:,:])
    z2[N1-k-1,:,:] = griddata(lonu.ravel(),latu.ravel(),z1.ravel(),vlon2,vlat2)
    Zlev2[N1-k-1]  = Zlev[k]

print ' \n' + '==> ' + '  INTERPOLATING V-Vel FROM Z --> S COORD ...\n' + ' '
V   = ztosigma(z2,ZsigV,Zlev2);

##############################################################################


### SOME PLOTTING ?? #########################################################

##############################################################################

# WRITING THE NETCDF FILE ####################################################
# Based on "ini_hydro.cdl" NETCDF sample structure

# some computings regarding netcdf variables:

N = int(run.klevels)
theta_s = run.theta_s
theta_b = run.theta_b
Mp, Lp = h2.shape
L  = Lp - 1
M  = Mp - 1
Np = N + 1

if run.spherical == 1:
    spherical = 'T'
else: 
    spherical = 'F'

ds       =  1.0 / N
lev      =  np.arange(1, N+1, 1)
sc       =  -1 + (lev-0.5)*ds
Ptheta   =  np.sinh(theta_s*sc) / np.sinh(theta_s)
Rtheta   =  np.tanh( theta_s*(sc+0.5) ) / ( 2* np.tanh(0.5*theta_s) ) - 0.5
Cs       =  (1-theta_b)*Ptheta + theta_b * Rtheta
scw      =  np.arange(-1, 0+ds, ds)
Pthetaw  =  np.sinh( theta_s*scw ) / np.sinh(theta_s)
Rthetaw  =  np.tanh( theta_s*(scw+0.5) ) / (2*np.tanh(0.5*theta_s)) - 0.5
Csw      =  (1-theta_b)*Pthetaw + theta_b*Rthetaw

print ' \n' + '==> ' + '  WRITING NETCDF INITIAL CONDITIONS FILE ...\n' + ' '

ncfile = netCDF4.Dataset(run.rundir + run.run_name + filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
ncfile.createDimension('s_rho', size=N)
ncfile.createDimension('s_w', size=Np)
ncfile.createDimension('tracer', size=2)
ncfile.createDimension('time', size=None)
ncfile.createDimension('one', size=1)

# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', filetypestr)
setattr(ncfile, 'title', run.ini_info)
setattr(ncfile, 'out_file', run.run_name + filenamestr)
setattr(ncfile, 'grd_file', run.run_name + '_grd.nc')
now = dt.datetime.now()
setattr(ncfile,'history',np.str(now))


# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'grid type logical switch')
setattr(ncfile.variables['spherical'], 'flag_values', 'T, F')
setattr(ncfile.variables['spherical'], 'flag_meanings', 'spherical, cartesian')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('Vtransform', 'd', dimensions=('one'))
setattr(ncfile.variables['Vtransform'], 'long_name',
    'vertical terrain-following transformation equation')
ncfile.variables['Vtransform'][:]  = run.vtransform

# ---------------------------------------------------------------------------
ncfile.createVariable('Vstretching', 'd', dimensions=('one'))
setattr(ncfile.variables['Vstretching'], 'long_name',
    'vertical terrain-following stretching function')
ncfile.variables['Vstretching'][:]  = run.vstretching

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_s', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_s'], 'long_name',
    'S-coordinate surface control parameter')
ncfile.variables['theta_s'][:] = run.theta_s

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_b', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_b'], 'long_name',
    'S-coordinate bottom control parameter')
ncfile.variables['theta_b'][:] = run.theta_b

# ---------------------------------------------------------------------------
ncfile.createVariable('Tcline', 'd', dimensions=('one'))
setattr(ncfile.variables['Tcline'], 'long_name',
    'S-coordinate surface/bottom layer width')
setattr(ncfile.variables['Tcline'], 'units', 'meter')
ncfile.variables['Tcline'][:] = run.tcline

# ---------------------------------------------------------------------------
ncfile.createVariable('hc', 'd', dimensions=('one'))
setattr(ncfile.variables['hc'],'long_name',
    'S-coordinate parameter, critical depth')
setattr(ncfile.variables['hc'], 'units', 'meter')
ncfile.variables['hc'][:] = run.hc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_rho', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['s_rho'], 'long_name', 'S-coordinate at RHO-points')
setattr(ncfile.variables['s_rho'], 'valid_min', -1.0)
setattr(ncfile.variables['s_rho'], 'valid_max', 0.0)
setattr(ncfile.variables['s_rho'], 'positive', 'up')
setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_rho'], 'formula_terms', 
    's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
ncfile.variables['s_rho'][:] = sc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['s_w'], 'long_name', 'S-coordinate at W-points')
setattr(ncfile.variables['s_w'], 'valid_min', -1.0)
setattr(ncfile.variables['s_w'], 'valid_max', 0.0)
setattr(ncfile.variables['s_w'], 'positive', 'up')
setattr(ncfile.variables['s_w'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_w'], 'formula_terms', 
    's: s_rho C: Cs_w eta: zeta depth: h depth_c: hc')
ncfile.variables['s_w'][:] = scw


# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_r', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['Cs_r'], 'long_name',
    'S-coordinate stretching curves at RHO-points')
setattr(ncfile.variables['Cs_r'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_r'], 'valid_max', 0.0)
ncfile.variables['Cs_r'][:] = Cs

# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['Cs_w'], 'long_name',
    'S-coordinate stretching curves at W-points')
setattr(ncfile.variables['Cs_w'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_w'], 'valid_max', 0.0)
ncfile.variables['Cs_w'][:] = Csw

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:] = h2

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree_east')
setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
ncfile.variables['lon_rho'][:] = grdfile.variables['lon_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree_north')
setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
ncfile.variables['lat_rho'][:] = grdfile.variables['lat_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree_east')
setattr(ncfile.variables['lon_u'], 'standard_name', 'longitude')
ncfile.variables['lon_u'][:] = grdfile.variables['lon_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree_north')
setattr(ncfile.variables['lat_u'], 'standard_name', 'latitude')
ncfile.variables['lat_u'][:] = grdfile.variables['lat_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree_east')
setattr(ncfile.variables['lon_v'], 'standard_name', 'lonitude')
ncfile.variables['lon_v'][:] = grdfile.variables['lon_v'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree_north')
setattr(ncfile.variables['lat_v'], 'standard_name', 'latitude')
ncfile.variables['lat_v'][:] = grdfile.variables['lat_v'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('ocean_time', 'd', dimensions=('time'))
setattr(ncfile.variables['ocean_time'], 'long_name',
    'time since initialization')
setattr(ncfile.variables['ocean_time'], 'units', 'seconds')
setattr(ncfile.variables['ocean_time'], 'calendar', '360.0 days / year')
ncfile.variables['ocean_time'] = run.time

# ---------------------------------------------------------------------------
ncfile.createVariable('zeta', 'd', dimensions=('time', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['zeta'], 'long_name', 'free-surface')
setattr(ncfile.variables['zeta'], 'units', 'meter')
setattr(ncfile.variables['zeta'], 'time', 'ocean_time')
setattr(ncfile.variables['zeta'], 'coordinates', 'lon_rho lat_rho ocean_time')
ncfile.variables['zeta'][:] = ZETA

# ---------------------------------------------------------------------------
ncfile.createVariable('ubar', 'd', dimensions=('time', 'eta_u', 'xi_u'))
setattr(ncfile.variables['ubar'], 'long_name', 
    'vertically integrated u-momentum component')
setattr(ncfile.variables['ubar'], 'units', 'meter second-1')
setattr(ncfile.variables['ubar'], 'time', 'ocean_time')
setattr(ncfile.variables['ubar'], 'coordinates', 'lon_u lat_u ocean_time')
ncfile.variables['ubar'][:] = UBAR

# ---------------------------------------------------------------------------
ncfile.createVariable('vbar', 'd', dimensions=('time', 'eta_v', 'xi_v'))
setattr(ncfile.variables['vbar'], 'long_name',
    'vertically integrated v-momentum component')
setattr(ncfile.variables['vbar'], 'units', 'meter second-1')
setattr(ncfile.variables['vbar'], 'time', 'ocean_time')
setattr(ncfile.variables['vbar'], 'coordinates', 'lon_v lat_v ocean_time')
ncfile.variables['vbar'][:] = VBAR

# ---------------------------------------------------------------------------
ncfile.createVariable('u', 'd', dimensions=('time', 's_rho', 'eta_u', 'xi_u'))
setattr(ncfile.variables['u'], 'long_name', 'u-momentum component')
setattr(ncfile.variables['u'], 'units', 'meter second-1')
setattr(ncfile.variables['u'], 'time', 'ocean_time')
setattr(ncfile.variables['u'], 'coordinates', 'lon_u lat_u s_rho ocean_time')
ncfile.variables['u'][:] = U

# ---------------------------------------------------------------------------
ncfile.createVariable('v', 'd', dimensions=('time', 's_rho', 'eta_v', 'xi_v'))
setattr(ncfile.variables['v'], 'long_name', 'v-momentum component')
setattr(ncfile.variables['v'], 'units', 'meter second-1')
setattr(ncfile.variables['v'], 'time', 'ocean_time')
setattr(ncfile.variables['v'], 'coordinates', 'lon_v lat_v s_rho ocean_time')
ncfile.variables['v'][:] = V

# ---------------------------------------------------------------------------
ncfile.createVariable('temp', 'd', dimensions=('time', 's_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['temp'], 'long_name', 'potential temperature')
setattr(ncfile.variables['temp'], 'units', 'celcius')
setattr(ncfile.variables['temp'], 'time', 'ocean_time')
setattr(ncfile.variables['temp'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['temp'][:] = TEMP

# ---------------------------------------------------------------------------
ncfile.createVariable('salt', 'd', dimensions=('time', 's_rho', 'eta_rho', 'xi_rho'))
setattr(ncfile.variables['salt'], 'long_name', 'salinity')
setattr(ncfile.variables['salt'], 'time', 'ocean_time')
setattr(ncfile.variables['salt'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
ncfile.variables['salt'][:] = SALT

ncfile.sync()

print ' \n' + '==> ' + '  ############################################  ...\n' + ' '
print ' \n' + '==> ' + '  INITIAL CONDITIONS FILE SUCCESSFULLY CREATED  ...\n' + ' '
print ' \n' + '==> ' + '  ############################################  ...\n' + ' '


### SOME PLOTTING  ###########################################################

print ' \n' + '==> ' + '  PLOTTING FIELDS FOR USER VERIFICATION  ...\n' + ' '

ncfile = netCDF4.Dataset(run.rundir + run.run_name + filenamestr, mode='r')
lon    = ncfile.variables['lon_rho'][:]
lat    = ncfile.variables['lat_rho'][:]
temp   = ncfile.variables['temp'][0, (N-1, N/2, 0), :, :]
salt   = ncfile.variables['salt'][0, (N-1, N/2, 0), :, :]

m = Basemap(projection='merc', llcrnrlat=run.latmin, urcrnrlat=run.latmax,
    llcrnrlon=run.lonmin, urcrnrlon=run.lonmax, lat_ts=0, resolution='l')

lon, lat = m(lon,lat) 

fig1 = plt.figure(1,figsize=(15,8),facecolor='w')

for k in (0, 1 ,2):
    plt.subplot(2,3,k+1)
    m.pcolormesh(lon, lat, np.squeeze(temp[k,:,:]))
    m.drawcoastlines(zorder=5)
    m.drawcountries(zorder=4)
    m.fillcontinents(color=('0.7'),lake_color='aqua',zorder=3)
    m.drawparallels(np.arange(run.latmin, run.latmax, 5),
	labels=[0, 0, 0, 0], dashes=[1,3], zorder=6)
    m.drawmeridians(np.arange(run.lonmin, run.lonmax, 5),
	labels=[0, 0, 0, 0], dashes=[1,3], zorder=7)

for k in (3, 4 ,5):
    plt.subplot(2,3,k+1)
    m.pcolormesh(lon, lat, np.squeeze(salt[k-3,:,:]))
    m.drawcoastlines(zorder=5)
    m.drawcountries(zorder=4)
    m.fillcontinents(color=('0.7'),lake_color='aqua',zorder=3)
    m.drawparallels(np.arange(run.latmin, run.latmax, 5),
	labels=[0, 0, 0, 0], dashes=[1,3], zorder=6)
    m.drawmeridians(np.arange(run.lonmin, run.lonmax, 5),
	labels=[0, 0, 0, 0], dashes=[1,3], zorder=7)


# print ' \n' + '==> ' + '  SAVING PNG FIGURE  ...\n' + ' '
# plt.savefig('figures/' + run.run_name + filenamestr[:-3] + '.png')

plt.show()

###############################################################################
del ncfile
filenamestr = '_bry.nc'
filetypestr = 'ROMS Boundary Conditions file'

# NORTH #####
# getting the northern slice to use as boundary condition
temp_north = TEMP[:,-1,:]; temp_north.shape = (1, run.klevels, Lp)
salt_north = SALT[:,-1,:]; salt_north.shape = (1, run.klevels, Lp)
zeta_north = ZETA[:,-1,:]; zeta_north.shape = (1, Lp)
u_north    =    U[:,-1,:];    u_north.shape = (1, run.klevels, L)
v_north    =    V[:,-1,:];    v_north.shape = (1, run.klevels, Lp)
# repeating as many times as the model will run
temp_north = temp_north.repeat(t.size, axis=0)
salt_north = salt_north.repeat(t.size, axis=0)
zeta_north = zeta_north.repeat(t.size, axis=0)
u_north    =    u_north.repeat(t.size, axis=0)
v_north    =    v_north.repeat(t.size, axis=0)

# EAST #######
# getting the eastern slice to use as boundary condition
temp_east = TEMP[:,:,-1]; temp_east.shape = (1, run.klevels, Mp)
salt_east = SALT[:,:,-1]; salt_east.shape = (1, run.klevels, Mp)
zeta_east = ZETA[:,:,-1]; zeta_east.shape = (1, Mp)
u_east    =    U[:,:,-1];    u_east.shape = (1, run.klevels, Mp)
v_east    =    V[:,:,-1];    v_east.shape = (1, run.klevels, M)
# repeating as many times as the model will run
temp_east = temp_east.repeat(t.size, axis=0)
salt_east = salt_east.repeat(t.size, axis=0)
zeta_east = zeta_east.repeat(t.size, axis=0)
u_east    =    u_east.repeat(t.size, axis=0)
v_east    =    v_east.repeat(t.size, axis=0)

# SOUTH #####
# getting the southern slice to use as boundary condition
temp_south = TEMP[:,1,:]; temp_south.shape = (1, run.klevels, Lp)
salt_south = SALT[:,1,:]; salt_south.shape = (1, run.klevels, Lp)
zeta_south = ZETA[:,1,:]; zeta_south.shape = (1, Lp)
u_south    =    U[:,1,:];    u_south.shape = (1, run.klevels, L)
v_south    =    V[:,1,:];    v_south.shape = (1, run.klevels, Lp)
# repeating as many times as the model will run
temp_south = temp_south.repeat(t.size, axis=0)
salt_south = salt_south.repeat(t.size, axis=0)
zeta_south = zeta_south.repeat(t.size, axis=0)
u_south    =    u_south.repeat(t.size, axis=0)
v_south    =    v_south.repeat(t.size, axis=0)


print ' \n' + '==> ' + '  WRITING NETCDF BOUNDARY CONDITIONS FILE ...\n' + ' '

ncfile = netCDF4.Dataset(run.rundir + run.run_name + filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
ncfile.createDimension('s_rho', size=N)
ncfile.createDimension('s_w', size=Np)
ncfile.createDimension('zeta_time', size=run.time)
ncfile.createDimension('v2d_time', size=run.time)
ncfile.createDimension('v3d_time', size=run.time)
ncfile.createDimension('temp_time', size=run.time)
ncfile.createDimension('salt_time', size=run.time)
ncfile.createDimension('one', size=1)

# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', filetypestr)
setattr(ncfile, 'title', run.ini_info)
setattr(ncfile, 'out_file', run.run_name + filenamestr)
setattr(ncfile, 'grd_file', run.run_name + '_grd.nc')
now = dt.datetime.now()
setattr(ncfile,'history',np.str(now))

# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'grid type logical switch')
setattr(ncfile.variables['spherical'], 'flag_values', 'T, F')
setattr(ncfile.variables['spherical'], 'flag_meanings', 'spherical, cartesian')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('Vtransform', 'd', dimensions=('one'))
setattr(ncfile.variables['Vtransform'], 'long_name',
    'vertical terrain-following transformation equation')
ncfile.variables['Vtransform'][:]  = run.vtransform

# ---------------------------------------------------------------------------
ncfile.createVariable('Vstretching', 'd', dimensions=('one'))
setattr(ncfile.variables['Vstretching'], 'long_name',
    'vertical terrain-following stretching function')
ncfile.variables['Vstretching'][:]  = run.vstretching

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_s', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_s'], 'long_name',
    'S-coordinate surface control parameter')
ncfile.variables['theta_s'][:] = run.theta_s

# ---------------------------------------------------------------------------
ncfile.createVariable('theta_b', 'd', dimensions=('one'))
setattr(ncfile.variables['theta_b'], 'long_name',
    'S-coordinate bottom control parameter')
ncfile.variables['theta_b'][:] = run.theta_b

# ---------------------------------------------------------------------------
ncfile.createVariable('Tcline', 'd', dimensions=('one'))
setattr(ncfile.variables['Tcline'], 'long_name',
    'S-coordinate surface/bottom layer width')
setattr(ncfile.variables['Tcline'], 'units', 'meter')
ncfile.variables['Tcline'][:] = run.tcline

# ---------------------------------------------------------------------------
ncfile.createVariable('hc', 'd', dimensions=('one'))
setattr(ncfile.variables['hc'],'long_name',
    'S-coordinate parameter, critical depth')
setattr(ncfile.variables['hc'], 'units', 'meter')
ncfile.variables['hc'][:] = run.hc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_rho', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['s_rho'], 'long_name', 'S-coordinate at RHO-points')
setattr(ncfile.variables['s_rho'], 'valid_min', -1.0)
setattr(ncfile.variables['s_rho'], 'valid_max', 0.0)
setattr(ncfile.variables['s_rho'], 'positive', 'up')
setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_rho'], 'formula_terms', 
    's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
ncfile.variables['s_rho'][:] = sc

# ---------------------------------------------------------------------------
ncfile.createVariable('s_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['s_w'], 'long_name', 'S-coordinate at W-points')
setattr(ncfile.variables['s_w'], 'valid_min', -1.0)
setattr(ncfile.variables['s_w'], 'valid_max', 0.0)
setattr(ncfile.variables['s_w'], 'positive', 'up')
setattr(ncfile.variables['s_w'], 'standard_name', 'ocean_s_coordinate_g1')
setattr(ncfile.variables['s_w'], 'formula_terms', 
    's: s_rho C: Cs_w eta: zeta depth: h depth_c: hc')
ncfile.variables['s_w'][:] = scw


# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_r', 'd', dimensions=('s_rho'))
setattr(ncfile.variables['Cs_r'], 'long_name',
    'S-coordinate stretching curves at RHO-points')
setattr(ncfile.variables['Cs_r'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_r'], 'valid_max', 0.0)
ncfile.variables['Cs_r'][:] = Cs

# ---------------------------------------------------------------------------
ncfile.createVariable('Cs_w', 'd', dimensions=('s_w'))
setattr(ncfile.variables['Cs_w'], 'long_name',
    'S-coordinate stretching curves at W-points')
setattr(ncfile.variables['Cs_w'], 'valid_min', -1.0)
setattr(ncfile.variables['Cs_w'], 'valid_max', 0.0)
ncfile.variables['Cs_w'][:] = Csw

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:] = h2

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree_east')
setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
ncfile.variables['lon_rho'][:] = grdfile.variables['lon_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree_north')
setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
ncfile.variables['lat_rho'][:] = grdfile.variables['lat_rho'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree_east')
setattr(ncfile.variables['lon_u'], 'standard_name', 'longitude')
ncfile.variables['lon_u'][:] = grdfile.variables['lon_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree_north')
setattr(ncfile.variables['lat_u'], 'standard_name', 'latitude')
ncfile.variables['lat_u'][:] = grdfile.variables['lat_u'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree_east')
setattr(ncfile.variables['lon_v'], 'standard_name', 'lonitude')
ncfile.variables['lon_v'][:] = grdfile.variables['lon_v'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree_north')
setattr(ncfile.variables['lat_v'], 'standard_name', 'latitude')
ncfile.variables['lat_v'][:] = grdfile.variables['lat_v'][:]

# ---------------------------------------------------------------------------
ncfile.createVariable('zeta_time', 'd', dimensions=('zeta_time'))
setattr(ncfile.variables['zeta_time'], 'long_name', 'free-surface time')
setattr(ncfile.variables['zeta_time'], 'units', 'days since 0000-01-01 00:00:00')
ncfile.variables['zeta_time'][:] = t

# ---------------------------------------------------------------------------
ncfile.createVariable('v3d_time', 'd', dimensions=('v3d_time'))
setattr(ncfile.variables['v3d_time'], 'long_name', '3D momentum time')
setattr(ncfile.variables['v3d_time'], 'units', 'days since 0000-01-01 00:00:00')
ncfile.variables['v3d_time'][:] = t

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_time', 'd', dimensions=('temp_time'))
setattr(ncfile.variables['temp_time'], 'long_name', 'potential temperature time')
setattr(ncfile.variables['temp_time'], 'units', 'days since 0000-01-01 00:00:00')
ncfile.variables['temp_time'][:] = t

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_time', 'd', dimensions=('salt_time'))
setattr(ncfile.variables['salt_time'], 'long_name', 'salinity time')
setattr(ncfile.variables['salt_time'], 'units', 'days since 0000-01-01 00:00:00')
ncfile.variables['salt_time'][:] = t

# ---------------------------------------------------------------------------
#ncfile.createVariable('u_west', 'd', dimensions=('v3d_time', 's_rho', 'eta_u'))
#setattr(ncfile.variables['u_west'], 'long_name', '3D u-momentum western boundary condition')
#setattr(ncfile.variables['u_west'], 'units', 'meter second-1')
#setattr(ncfile.variables['u_west'], 'time', 'v3d_time')
#ncfile.variables['u_west'][:] = u_west

# ---------------------------------------------------------------------------
ncfile.createVariable('u_east', 'd', dimensions=('v3d_time', 's_rho', 'eta_u'))
setattr(ncfile.variables['u_east'], 'long_name', '3D u-momentum eastern boundary condition')
setattr(ncfile.variables['u_east'], 'units', 'meter second-1')
setattr(ncfile.variables['u_east'], 'time', 'v3d_time')
ncfile.variables['u_east'][:] = u_east

# ---------------------------------------------------------------------------
ncfile.createVariable('u_south', 'd', dimensions=('v3d_time', 's_rho', 'xi_u'))
setattr(ncfile.variables['u_south'], 'long_name', '3D u-momentum southern boundary condition')
setattr(ncfile.variables['u_south'], 'units', 'meter second-1')
setattr(ncfile.variables['u_south'], 'time', 'v3d_time')
ncfile.variables['u_south'][:] = u_south

# ---------------------------------------------------------------------------
ncfile.createVariable('u_north', 'd', dimensions=('v3d_time', 's_rho', 'xi_u'))
setattr(ncfile.variables['u_north'], 'long_name', '3D u-momentum northern boundary condition')
setattr(ncfile.variables['u_north'], 'units', 'meter second-1')
setattr(ncfile.variables['u_north'], 'time', 'v3d_time')
ncfile.variables['u_north'][:] = u_north

## ---------------------------------------------------------------------------
#ncfile.createVariable('v_west', 'd', dimensions=('v3d_time', 's_rho', 'eta_v'))
#setattr(ncfile.variables['v_west'], 'long_name', '3D v-momentum western boundary condition')
#setattr(ncfile.variables['v_west'], 'units', 'meter second-1')
#setattr(ncfile.variables['v_west'], 'time', 'v3d_time')
#ncfile.variables['v_west'][:] = v_west

# ---------------------------------------------------------------------------
ncfile.createVariable('v_east', 'd', dimensions=('v3d_time', 's_rho', 'eta_v'))
setattr(ncfile.variables['v_east'], 'long_name', '3D v-momentum eastern boundary condition')
setattr(ncfile.variables['v_east'], 'units', 'meter second-1')
setattr(ncfile.variables['v_east'], 'time', 'v3d_time')
ncfile.variables['v_east'][:] = v_east

# ---------------------------------------------------------------------------
ncfile.createVariable('v_south', 'd', dimensions=('v3d_time', 's_rho', 'xi_v'))
setattr(ncfile.variables['v_south'], 'long_name', '3D v-momentum sovthern boundary condition')
setattr(ncfile.variables['v_south'], 'units', 'meter second-1')
setattr(ncfile.variables['v_south'], 'time', 'v3d_time')
ncfile.variables['v_south'][:] = v_south

# ---------------------------------------------------------------------------
ncfile.createVariable('v_north', 'd', dimensions=('v3d_time', 's_rho', 'xi_v'))
setattr(ncfile.variables['v_north'], 'long_name', '3D v-momentum northern boundary condition')
setattr(ncfile.variables['v_north'], 'units', 'meter second-1')
setattr(ncfile.variables['v_north'], 'time', 'v3d_time')
ncfile.variables['v_north'][:] = v_north

# ---------------------------------------------------------------------------
#ncfile.createVariable('temp_west', 'd', dimensions=('temp_time', 's_rho', 'eta_rho'))
#setattr(ncfile.variables['temp_west'], 'long_name', 'potential temperature western boundary condition')
#setattr(ncfile.variables['temp_west'], 'units', 'celcius')
#setattr(ncfile.variables['temp_west'], 'time', 'temp_time')
#ncfile.variables['temp_west'][:] = temp_west

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_east', 'd', dimensions=('temp_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['temp_east'], 'long_name', 'potential temperature eastern boundary condition')
setattr(ncfile.variables['temp_east'], 'units', 'celcius')
setattr(ncfile.variables['temp_east'], 'time', 'temp_time')
ncfile.variables['temp_east'][:] = temp_east

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_south', 'd', dimensions=('temp_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['temp_south'], 'long_name', 'potential temperature southern boundary condition')
setattr(ncfile.variables['temp_south'], 'units', 'celcius')
setattr(ncfile.variables['temp_south'], 'time', 'temp_time')
ncfile.variables['temp_south'][:] = temp_south

# ---------------------------------------------------------------------------
ncfile.createVariable('temp_north', 'd', dimensions=('temp_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['temp_north'], 'long_name', 'potential temperature northern boundary condition')
setattr(ncfile.variables['temp_north'], 'units', 'celcius')
setattr(ncfile.variables['temp_north'], 'time', 'temp_time')
ncfile.variables['temp_north'][:] = temp_north

## ---------------------------------------------------------------------------
#ncfile.createVariable('salt_west', 'd', dimensions=('salt_time', 's_rho', 'eta_rho'))
#setattr(ncfile.variables['salt_west'], 'long_name', 'salinity western boundary condition')
#setattr(ncfile.variables['salt_west'], 'time', 'salt_time')
#ncfile.variables['salt_west'][:] = salt_west

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_east', 'd', dimensions=('salt_time', 's_rho', 'eta_rho'))
setattr(ncfile.variables['salt_east'], 'long_name', 'salinity eastern boundary condition')
setattr(ncfile.variables['salt_east'], 'time', 'salt_time')
ncfile.variables['salt_east'][:] = salt_east

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_south', 'd', dimensions=('salt_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['salt_south'], 'long_name', 'salinity southern boundary condition')
setattr(ncfile.variables['salt_south'], 'time', 'salt_time')
ncfile.variables['salt_south'][:] = salt_south

# ---------------------------------------------------------------------------
ncfile.createVariable('salt_north', 'd', dimensions=('salt_time', 's_rho', 'xi_rho'))
setattr(ncfile.variables['salt_north'], 'long_name', 'salinity northern boundary condition')
setattr(ncfile.variables['salt_north'], 'time', 'salt_time')
ncfile.variables['salt_north'][:] = salt_north

# ---------------------------------------------------------------------------
#ncfile.createVariable('zeta_west', 'd', dimensions=('zeta_time', 's_rho', 'eta_rho'))
#setattr(ncfile.variables['zeta_west'], 'long_name', 'free-surface western boundary condition')
#setattr(ncfile.variables['zeta_west'], 'units', 'meter')
#setattr(ncfile.variables['zeta_west'], 'time', 'zeta_time')
#ncfile.variables['zeta_west'][:] = zeta_west

# ---------------------------------------------------------------------------
ncfile.createVariable('zeta_east', 'd', dimensions=('zeta_time', 'eta_rho'))
setattr(ncfile.variables['zeta_east'], 'long_name', 'free-surface eastern boundary condition')
setattr(ncfile.variables['zeta_east'], 'units', 'meter')
setattr(ncfile.variables['zeta_east'], 'time', 'zeta_time')
ncfile.variables['zeta_east'][:] = zeta_east

# ---------------------------------------------------------------------------
ncfile.createVariable('zeta_south', 'd', dimensions=('zeta_time', 'xi_rho'))
setattr(ncfile.variables['zeta_south'], 'long_name', 'free-surface southern boundary condition')
setattr(ncfile.variables['zeta_south'], 'units', 'meter')
setattr(ncfile.variables['zeta_south'], 'time', 'zeta_time')
ncfile.variables['zeta_south'][:] = zeta_south

# ---------------------------------------------------------------------------
ncfile.createVariable('zeta_north', 'd', dimensions=('zeta_time', 'xi_rho'))
setattr(ncfile.variables['zeta_north'], 'long_name', 'free-surface northern boundary condition')
setattr(ncfile.variables['zeta_north'], 'units', 'meter')
setattr(ncfile.variables['zeta_north'], 'time', 'zeta_time')
ncfile.variables['zeta_north'][:] = zeta_north

ncfile.sync()

print ' \n' + '==> ' + '  #############################################  ...\n' + ' '
print ' \n' + '==> ' + '  BOUNDARY CONDITIONS FILE SUCCESSFULLY CREATED  ...\n' + ' '
print ' \n' + '==> ' + '  #############################################  ...\n' + ' '

###########

print ' \n' + '==> ' + '  ############################################  ...\n' + ' '
print ' \n' + '==> ' + '                  !! DONE !!                    ...\n' + ' '
print ' \n' + '==> ' + '  ############################################  ...\n' + ' '

##############################################################################

