#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script to download HYCOM outputs via opendap and seve as netcdf
#  HYCOM + NCODA Global 1/12 Analysis (expt_60.5) (assimilates data) Nov-2003 to Dec-2006
#  url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5';
#  
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Apr, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from matplotlib.mlab import find
import datetime as dt
import netCDF4 as nc
from pydap.client import open_url

def near(x,x0):
    """
    Find the index where x has the closer value to x0
    """
    
    dx = x - x0
    dx = np.abs(dx)
    fn = np.where( dx == dx.min() )
    fn = fn[0][0]
    
    return fn
    
###########################################################################

print 'Accessing dods server'
dataset = open_url('http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5')

print 'Downloading lon, lat, depth, time arrays'
tt  = dataset.MT[:]
depth  = dataset['Depth'][:]
lon = dataset.Longitude.array[:]; lon = lon-360
lat = dataset.Latitude.array[:]

# obs: hycom dates are referenced to year 1900
tstart = dt.datetime(2003-1900, 1, 1); tstart = date2num(tstart)
tend   = dt.datetime(2006-1900, 1, 1); tend   = date2num(tend)
fk  = near(depth, 60) 
fj  = np.where( (lon[1,:] >= -60) & (lon[1,:] <= 30) )
fi  = np.where( (lat[:,1] >= -50) & (lat[:,1] <= -5) )
fj1 = fj[0][0]; fj2 = fj[0][1]
fi1 = fi[0][0]; fi2 = fi[0][1]

ft  = find( (tt > tstart) & (tt < tend) )  
tt2 = tt[ft]
lon = dataset.Longitude.array[ fi1:fi2, fj1:fj2 ]; lon = lon-360
lat = dataset.Latitude.array[ fi1:fi2, fj1:fj2 ]
depth  = dataset['Depth'][fk]

# downloading data from the server according to the values prescribed before:
uhyc = dataset['u']
vhyc = dataset['v']
#thyc = dataset['temperature']
#shyc = dataset['salinity']
hhyc = dataset['ssh']

print 'Downloading U, V, T, S arrays for each time-step and storing in netCDF file:'
# downloading data from the server according to the values prescribed before

print '    Initializing netCDF file to record data'
ncfile = nc.Dataset('outputs/hycom_glba0.08_exp60.5_60m.nc',
            mode='w', clobber='true', format='NETCDF3_CLASSIC')

# dimensions
ncfile.createDimension('one', size=1)
ncfile.createDimension('time', size=tt2.size)
ncfile.createDimension('Y', size = fi2 - fi1)
ncfile.createDimension('X', size = fj2 - fj1)

# global attributes  
setattr(ncfile, 'title','HYCOM GLBa0.08 expt 60.5, 2000-2006 daily fields at 60 m')
now = dt.datetime.now()
setattr(ncfile,'created',np.str(now))

# variables
ncfile.createVariable('u', 'd', dimensions=('time', 'Y', 'X'))
setattr(ncfile.variables['u'], 'long_name', 'U-velocity [m/s]')

ncfile.createVariable('v', 'd', dimensions=('time', 'Y', 'X'))
setattr(ncfile.variables['v'], 'long_name', 'V-velocity [m/s]')

#ncfile.createVariable('temp', 'd', dimensions=('time', 'Y', 'X'))
#setattr(ncfile.variables['temp'], 'long_name', 'Temperature [degree C]')

#ncfile.createVariable('salt', 'd', dimensions=('time', 'Y', 'X'))
#setattr(ncfile.variables['salt'], 'long_name', 'Salinity [ ]')

ncfile.createVariable('ssh', 'd', dimensions=('time', 'Y', 'X'))
setattr(ncfile.variables['ssh'], 'long_name', 'Sea Surface Height [m]')

ncfile.createVariable('lon', 'd', dimensions=('Y', 'X'))
setattr(ncfile.variables['lon'], 'long_name', 'Longitude')
ncfile.variables['lon'][:] = lon

ncfile.createVariable('lat', 'd', dimensions=('Y', 'X'))
setattr(ncfile.variables['lat'], 'long_name', 'Latitude')
ncfile.variables['lat'][:] = lat    
    
c = 0
for k in tt2:
    ft  = find(tt == k)
    time = num2date(k + 693961)
    print '    HYCOM time step: ', time

    U   = uhyc.array[ ft, fk, fi1:fi2, fj1:fj2 ] 
    ncfile.variables['u'][c,...] = U

    V   = vhyc.array[ ft, fk, fi1:fi2, fj1:fj2 ]
    ncfile.variables['v'][c,...] = V

    #T   = thyc.array[ ft, fk, fi1:fi2, fj1:fj2 ]
    #ncfile.variables['temp'][c,...] = T

    #S   = shyc.array[ ft, fk, fi1:fi2, fj1:fj2 ]
    #ncfile.variables['salt'][c,...] = S
    
    H   = hhyc.array[ ft, fk, fi1:fi2, fj1:fj2 ]
    ncfile.variables['ssh'][c,...] = H
	
    del U, V, H, ncfile
    c   = c + 1

ncfile.sync()   
 
