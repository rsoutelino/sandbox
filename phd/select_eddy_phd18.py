#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  plot_roms.py -- Visualization of ROMS fields
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Jun, 2013
###################################################################################
import sys
import numpy as np
import scipy.io as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near
from romslab import subset



def near(x,x0):
    """
    Find the index where x has the closer value to x0
    """
    
    dx = x - x0
    dx = abs(dx)
    fn = np.where( dx == dx.min() )
    
    return fn[0][0]


# PYTHON SCRIPT START #####################################################

# Input parameters
script, expt, zlev, eddy = sys.argv
zlev = int(zlev)

# expt = 'phd15'
# zlev = 0
# eddy = 'AbE'

zlev = int(zlev)


if eddy == "IE":
    lims = [-38.53, -36.59, -15.689, -13.23]  # limits for Ilheus Eddy
elif eddy == "RCE":
    lims = [-38.34, -36.53, -17.66, -16.16]      # limits for Royal-Charlotte eddy
elif eddy == "AbE":
    lims = [-37.65, -35.78, -20.51, -17.91]   # limits for Abrolhos Eddy


print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
roms = run_setup(expt + '_run.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grd  = netCDF4.Dataset(roms.rundir + roms.run_name + '_grd.nc')

# assigning some variables from grid file
lonr  = grd.variables['lon_rho'][:]
latr  = grd.variables['lat_rho'][:]
latu  = grd.variables['lat_u'][:]
latv  = grd.variables['lat_v'][:]
lonu  = grd.variables['lon_u'][:]
lonv  = grd.variables['lon_v'][:]
h     = grd.variables['h'][:]

print ' \n' + '==> ' + '  READING CHOSEN ROMS OUTPUT NETCDF FILE ...\n' + ' '
avg   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg.nc')

lm, km, im, jm = avg.variables['temp'].shape

tmp = avg.variables['temp'][0,0,...]

lon, lat, tmp = subset(lonr, latr, tmp, lims[0], lims[1], lims[2], lims[3])


# AVG  ===============================================================

zu = get_depths(avg, grd, 0, 'u')
zv = get_depths(avg, grd, 0, 'v')
zt = get_depths(avg, grd, 0, 't')

times = np.arange(1800, lm, 10)
U = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
V = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
T = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)


print "\n\n      Slicing %s eddy, at %s m, in AVG.... may take a while.......\n\n" %(eddy, zlev)

if zlev == 0:
    count = 0
    for l in times:
        print "AVG - step %s of %s, at %s m" %(l, times[-1], zlev)
        u = avg.variables['u'][l,-1,...]
        v = avg.variables['v'][l,-1,...]
        t = avg.variables['temp'][l,-1,...]
        x, y, u = subset(lonu, latu, u, lims[0], lims[1], lims[2], lims[3])
        x, y, v = subset(lonu, latu, v, lims[0], lims[1], lims[2], lims[3])
        x, y, t = subset(lonu, latu, t, lims[0], lims[1], lims[2], lims[3])
        U[count,...] = u
        V[count,...] = v
        T[count,...] = t
        count += 1
else:
    count = 0
    for l in times:
        print "AVG - step %s of %s: S --> Z, at %s m" %(l, times[-1], zlev)
        u = avg.variables['u'][l,:,...]
        v = avg.variables['v'][l,:,...]
        t = avg.variables['temp'][l,:,...]
        utmp = lonu*0
        vtmp = lonv*0
        ttmp = lonr*0

        for i in range(0, im):
            for j in range(0, jm-1):
                utmp[i, j] = np.interp(-zlev, zu[:, i, j], u[:, i, j] )

        for i in range(0, im-1):
            for j in range(0, jm):
                vtmp[i, j] = np.interp(-zlev, zv[:, i, j], v[:, i, j] )

        for i in range(0, im):
            for j in range(0, jm):
                ttmp[i, j] = np.interp(-zlev, zt[:, i, j], t[:, i, j] )
        
        x, y, utmp = subset(lonu, latu, utmp, lims[0], lims[1], lims[2], lims[3])
        x, y, vtmp = subset(lonu, latu, vtmp, lims[0], lims[1], lims[2], lims[3])
        x, y, ttmp = subset(lonu, latu, ttmp, lims[0], lims[1], lims[2], lims[3])
        U[count,...] = utmp
        V[count,...] = vtmp
        T[count,...] = ttmp
        count += 1
        


U = np.ma.masked_where(np.abs(U) > 2, U)
V = np.ma.masked_where(np.abs(V) > 2, V)
T = np.ma.masked_where(np.abs(T) > 40, T)



mdict = {'lon':lon, 'lat':lat, 'u':U, 'v':V, 'temp':T}
filename = "eddies_subsets/%s_%s_subset_%sm.mat" %(expt, eddy, zlev)
sp.savemat(filename, mdict)


m = Basemap(projection='cyl', llcrnrlat = latr.min(), urcrnrlat = latr.max(), 
            llcrnrlon = lonr.min(), urcrnrlon = lonr.max(),
            lat_ts=0, resolution='i')

# plt.figure(figsize=(10,10))
# m.pcolormesh(lonr, latr, h, cmap=plt.cm.Blues)
# m.quiver(lon[::2,::2], lat[::2,::2], U.mean(axis=0)[::2,::2], V.mean(axis=0)[::2,::2])
# m.fillcontinents()
# plt.title("%s Eddy subset on level: %s m" %(eddy,zlev) )
# plt.show()

