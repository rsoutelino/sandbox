#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  plot_roms.py -- Visualization of ROMS fields
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Jun, 2013
###################################################################################

import numpy as np
import scipy.io as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near



def near(x,x0):
    """
    Find the index where x has the closer value to x0
    """
    
    dx = x - x0
    dx = abs(dx)
    fn = np.where( dx == dx.min() )
    
    return fn[0][0]


# PYTHON SCRIPT START ######################################################

# Input parameters
expt = 'phd16'
zlev = 2000


print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
roms = run_setup(expt + '_run.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grd  = netCDF4.Dataset(roms.rundir + roms.run_name + '_grd.nc')

# assigning some variables from grid file
lon   = grd.variables['lon_rho'][:]
lat   = grd.variables['lat_rho'][:]
latu  = grd.variables['lat_u'][:]
latv  = grd.variables['lat_v'][:]
lonu  = grd.variables['lon_u'][:]
lonv  = grd.variables['lon_v'][:]
h     = grd.variables['h'][:]

print ' \n' + '==> ' + '  READING CHOSEN ROMS OUTPUT NETCDF FILE ...\n' + ' '
avg1   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-1.nc')
avg2   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-2.nc')
avg3   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-3.nc')



lm1, km, im, jm = avg1.variables['temp'].shape
lm2, km, im, jm = avg2.variables['temp'].shape
lm3, km, im, jm = avg3.variables['temp'].shape


print "\n\n      Computing run average .... may take a while.......\n\n"

if zlev == 0:
    print "Computing average for AVG 1\n\n" 
    umean1 = avg1.variables['u'][1800::10,-1,...].mean(axis=0)
    vmean1 = avg1.variables['v'][1800::10,-1,...].mean(axis=0)

    print "Computing average for AVG 2\n\n" 
    umean2 = avg2.variables['u'][::10,-1,...].mean(axis=0)
    vmean2 = avg2.variables['v'][::10,-1,...].mean(axis=0)

    print "Computing average for AVG 3\n\n" 
    umean3 = avg3.variables['u'][::10,-1,...].mean(axis=0)
    vmean3 = avg3.variables['v'][::10,-1,...].mean(axis=0)  

    umean = (umean1 + umean2 + umean3) / 3
    vmean = (vmean1 + vmean2 + vmean3) / 3

    umean = griddata(lonu.ravel(), latu.ravel(), umean.ravel(), lon, lat)
    vmean = griddata(lonv.ravel(), latv.ravel(), vmean.ravel(), lon, lat) 

else:
    zu = get_depths(avg1, grd, 0, 'u')
    zv = get_depths(avg1, grd, 0, 'v')
    umean1 = np.zeros( (km, im, jm-1) )
    vmean1 = np.zeros( (km, im-1, jm) )
    umean2 = np.zeros( (km, im, jm-1) )
    vmean2 = np.zeros( (km, im-1, jm) )
    umean3 = np.zeros( (km, im, jm-1) )
    vmean3 = np.zeros( (km, im-1, jm) )
    umean  = np.zeros( (im, jm-1) )
    vmean  = np.zeros( (im-1, jm) )
    count = 0

    print "\nComputing average for AVG 1" 
    for l in np.arange(1800, lm1, 10):
        print "    Day %s of %s" %(l, lm1)
        count += 1
        umean1 += avg1.variables['u'][l,...]
        vmean1 += avg1.variables['v'][l,...]

    umean1 = umean1 / count
    vmean1 = vmean1 / count


    count = 0

    print "\n\nComputing average for AVG 2"
    for l in np.arange(0, lm2, 10):
        print "    Day %s of %s" %(l, lm2)
        count += 1
        count += 1
        umean2 += avg2.variables['u'][l,...]
        vmean2 += avg2.variables['v'][l,...]
   
    umean2 = umean2 / count
    vmean2 = vmean2 / count


    count = 0

    print "\n\nComputing average for AVG 3"
    for l in np.arange(0, lm3, 10):
        print "    Day %s of %s" %(l, lm3)
        count += 1
        count += 1
        umean3 += avg3.variables['u'][l,...]
        vmean3 += avg3.variables['v'][l,...]
   
    umean3 = umean3 / count
    vmean3 = vmean3 / count


    u = (umean1 + umean2 + umean3) / 3
    v = (vmean1 + vmean2 + vmean3) / 3
    u = np.ma.masked_where(np.abs(u) > 2, u)
    v = np.ma.masked_where(np.abs(v) > 2, v)

    print "\n\n\n    S --> Z:"
    for i in range(0, im):
        for j in range(0, jm-1):
            umean[i,j] = np.interp(-zlev, zu[:, i, j], u[:, i, j] )

    for i in range(0, im-1):
        for j in range(0, jm):
            vmean[i,j] = np.interp(-zlev, zv[:, i, j], v[:, i, j] )


    umean = griddata(lonu.ravel(), latu.ravel(), umean.ravel(), lon, lat) 
    vmean = griddata(lonv.ravel(), latv.ravel(), vmean.ravel(), lon, lat) 
    umean = np.ma.masked_where(h < zlev, umean)
    vmean = np.ma.masked_where(h < zlev, vmean) 


umean = np.ma.masked_where(np.abs(umean) > 2, umean)
vmean = np.ma.masked_where(np.abs(vmean) > 2, vmean)

mdict = {'lon':lon,'lat':lat,'umean':umean,'vmean':vmean}
filename = "averages/%s_run-average-vels_%sm.mat" %(expt, zlev)
sp.savemat(filename, mdict)

if zlev != 0:
    mdict2 = {'lon':lon,'lat':lat,'umean':u,'vmean':v}
    filename2 = "averages/%s_run-average-vels.mat" %expt
    sp.savemat(filename2, mdict2)



m = Basemap(projection='cyl', llcrnrlat = lat.min(), urcrnrlat = lat.max(), 
            llcrnrlon = lon.min(), urcrnrlon = lon.max(),
            lat_ts=0, resolution='l')

plt.figure(figsize=(10,10))
m.quiver(lon[::4,::4], lat[::4,::4], umean[::4,::4], vmean[::4,::4])
m.fillcontinents()
plt.title("Run-averaged velocities: %s m" %zlev)
plt.show()

print "\n\n     Finished processing %s at %sm " %(expt, zlev)