#!/usr/bin/env python
#
# Creates forcing fields netCDF file for ROMS
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
from roms_setup import run_setup, wind_stress

#####################################################################

# SCRIPT START ######################################################

# Basic Settings:

filenamestr = '_frc.nc'
filetypestr = 'ROMS Forcing file'

# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# OA-created netcdf initial T, S file 
# grid netcdf file

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
run = run_setup('../phd2_run.setup')

print ' \n' + '==> ' + '  READING ECMWF NETCDF FILE ...\n' + ' '
datafile = sp.loadmat(run.datadir + run.frc_filename)

# assigning some variables from data file
lon = datafile.pop('lon')
lat = datafile.pop('lat')
u10 = datafile.pop('u')
v10 = datafile.pop('v')

print ' \n' + '==> ' + '  COMPUTING WIND STRESS ...\n' + ' '
sustr, svstr = wind_stress(u10, v10)

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(run.rundir + run.run_name + '_grd.nc')

# assigning some variables from grid file
rlon2  = grdfile.variables['lon_rho'][:]
rlat2  = grdfile.variables['lat_rho'][:]
lonu  = grdfile.variables['lon_u'][:]
latu  = grdfile.variables['lat_u'][:]
lonv  = grdfile.variables['lon_v'][:]
latv  = grdfile.variables['lat_v'][:]
h2    = grdfile.variables['h'][:]

# DOING COMPUTATIONS TO INTERPOLATE THE FIELDS TO ROMS GRID #########
t = np.arange(0, run.time);
Jrho, Irho = rlon2.shape
Mr2, Lr2   = rlon2.shape
Lu2 = Lr2-1; Mu2 = Mr2
Lv2 = Lr2;   Mv2 = Mr2-1

### Interpolating Wind Field to ROMS 2D grid ###############################
# Currently doing only TEMP and SALT. ZETA, UBAR, VBAR, U, V are zero arrays

print ' \n' + '==> ' + '  INTERPOLATING WIND ...\n' + ' '
sustr = griddata(lon.ravel(), lat.ravel(), sustr.ravel(), lonu, latu) 
svstr = griddata(lon.ravel(), lat.ravel(), svstr.ravel(), lonv, latv) 

# repeating as many times as the model will run
sustr.shape = (1, Mu2, Lu2)
svstr.shape = (1, Mv2, Lv2)
sustr = sustr.repeat(t.size, axis=0)
svstr = svstr.repeat(t.size, axis=0)

# WRITING THE NETCDF FILE ####################################################
# Based on "frc_uvstress.cdl" NETCDF sample structure

# some computings regarding netcdf variables:
Mp, Lp = h2.shape
L  = Lp - 1
M  = Mp - 1

print ' \n' + '==> ' + '  WRITING NETCDF FORCING FILE ...\n' + ' '

ncfile = netCDF4.Dataset(run.rundir + run.run_name + filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
ncfile.createDimension('sms_time', size=run.time)
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
ncfile.createVariable('sms_time', 'd', dimensions=('sms_time'))
setattr(ncfile.variables['sms_time'], 'long_name', 'surface momentum stress time')
setattr(ncfile.variables['sms_time'], 'units', 'days since 0000-01-01 00:00:00')
setattr(ncfile.variables['sms_time'], 'interval', '360.0 daily data')
setattr(ncfile.variables['sms_time'], 'start_date', '01-Jan-0000')
setattr(ncfile.variables['sms_time'], 'end_date', '10-Feb-0000')
ncfile.variables['sms_time'][:] = t

# ---------------------------------------------------------------------------
ncfile.createVariable('sustr', 'd', dimensions=('sms_time', 'eta_u', 'xi_u'))
setattr(ncfile.variables['sustr'], 'long_name', 'surface u-momentum stress')
setattr(ncfile.variables['sustr'], 'units', 'Newton meter-2')
setattr(ncfile.variables['sustr'], 'time', 'sms_time')
ncfile.variables['sustr'][:] = sustr

# ---------------------------------------------------------------------------
ncfile.createVariable('svstr', 'd', dimensions=('sms_time', 'eta_v', 'xi_v'))
setattr(ncfile.variables['svstr'], 'long_name', 'surface v-momentum stress')
setattr(ncfile.variables['svstr'], 'units', 'Newton meter-2')
setattr(ncfile.variables['svstr'], 'time', 'sms_time')
ncfile.variables['svstr'][:] = svstr

ncfile.sync()

print ' \n' + '==> ' + '  ################################# ...\n' + ' '
print ' \n' + '==> ' + '  FORCING FILE SUCCESSFULLY CREATED ...\n' + ' '
print ' \n' + '==> ' + '  ################################# ...\n' + ' '



