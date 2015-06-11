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
avg1   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-1.nc')
avg2   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-2.nc')
avg3   = netCDF4.Dataset(roms.rundir + roms.run_name + '_avg-3.nc')

lm1, km, im, jm = avg1.variables['temp'].shape
lm2, km, im, jm = avg2.variables['temp'].shape
lm3, km, im, jm = avg3.variables['temp'].shape



tmp = avg1.variables['temp'][0,0,...]

lon, lat, tmp = subset(lonr, latr, tmp, lims[0], lims[1], lims[2], lims[3])


# AVG1 output ===============================================================

zu = get_depths(avg1, grd, 0, 'u')
zv = get_depths(avg1, grd, 0, 'v')
zt = get_depths(avg1, grd, 0, 't')

times = np.arange(1800, lm1, 10)
U1 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
V1 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
T1 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)


print "\n\n      Slicing %s eddy, at %s m, in AVG-1.... may take a while.......\n\n" %(eddy, zlev)

if zlev == 0:
    count = 0
    for l in times:
        print "AVG-1 - step %s of %s, at %s m" %(l, times[-1], zlev)
        u = avg1.variables['u'][l,-1,...]
        v = avg1.variables['v'][l,-1,...]
        t = avg1.variables['temp'][l,-1,...]
        x, y, u = subset(lonu, latu, u, lims[0], lims[1], lims[2], lims[3])
        x, y, v = subset(lonu, latu, v, lims[0], lims[1], lims[2], lims[3])
        x, y, t = subset(lonu, latu, t, lims[0], lims[1], lims[2], lims[3])
        U1[count,...] = u
        V1[count,...] = v
        T1[count,...] = t
        count += 1
else:
    count = 0
    for l in times:
        print "AVG-1 - step %s of %s: S --> Z, at %s m" %(l, times[-1], zlev)
        u = avg1.variables['u'][l,:,...]
        v = avg1.variables['v'][l,:,...]
        t = avg1.variables['temp'][l,:,...]
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
        U1[count,...] = utmp
        V1[count,...] = vtmp
        T1[count,...] = ttmp
        count += 1
        

 # AVG2 output ===============================================================   


zu = get_depths(avg2, grd, 0, 'u')
zv = get_depths(avg2, grd, 0, 'v')
zt = get_depths(avg2, grd, 0, 't')

times = np.arange(0, lm2, 10)
U2 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
V2 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
T2 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)

print "\n\n      Slicing the %s eddy, at %s m, in AVG-2.... may take a while.......\n\n" %(eddy, zlev)

if zlev == 0:
    count = 0
    for l in times:
        print "AVG-2 - step %s of %s, at %s m" %(l, times[-1], zlev)
        u = avg2.variables['u'][l,-1,...]
        v = avg2.variables['v'][l,-1,...]
        t = avg2.variables['temp'][l,-1,...]
        x, y, u = subset(lonu, latu, u, lims[0], lims[1], lims[2], lims[3])
        x, y, v = subset(lonu, latu, v, lims[0], lims[1], lims[2], lims[3])
        x, y, t = subset(lonu, latu, t, lims[0], lims[1], lims[2], lims[3])
        U2[count,...] = u
        V2[count,...] = v
        T2[count,...] = t
        count += 1
else:
    count = 0
    for l in times:
        print "AVG-2 - step %s of %s: S --> Z, at %s m" %(l, times[-1], zlev)
        u = avg2.variables['u'][l,:,...]
        v = avg2.variables['v'][l,:,...]
        t = avg2.variables['temp'][l,:,...]
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
        U2[count,...] = utmp
        V2[count,...] = vtmp
        T2[count,...] = ttmp
        count += 1


# AVG3 output ===============================================================   


zu = get_depths(avg3, grd, 0, 'u')
zv = get_depths(avg3, grd, 0, 'v')
zt = get_depths(avg3, grd, 0, 't')
times = np.arange(0, lm3, 10)
U3 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
V3 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)
T3 = np.zeros((times.size, lon.shape[0], lon.shape[1]), dtype=np.float64)


print "\n\n      Slicing the %s eddy, at %s m, in AVG-2.... may take a while.......\n\n" %(eddy, zlev)

if zlev == 0:
    count = 0
    for l in times:
        print "AVG-3 - step %s of %s, at %s m" %(l, times[-1], zlev)
        u = avg3.variables['u'][l,-1,...]
        v = avg3.variables['v'][l,-1,...]
        t = avg3.variables['temp'][l,-1,...]
        x, y, u = subset(lonu, latu, u, lims[0], lims[1], lims[2], lims[3])
        x, y, v = subset(lonu, latu, v, lims[0], lims[1], lims[2], lims[3])
        x, y, t = subset(lonu, latu, t, lims[0], lims[1], lims[2], lims[3])
        U3[count,...] = u
        V3[count,...] = v
        T3[count,...] = t
        count += 1
else:
    count = 0
    for l in times:
        print "AVG-3 - step %s of %s: S --> Z, at %s m" %(l, times[-1], zlev)
        u = avg3.variables['u'][l,:,...]
        v = avg3.variables['v'][l,:,...]
        t = avg3.variables['temp'][l,:,...]
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
        U3[count,...] = utmp
        V3[count,...] = vtmp
        T3[count,...] = ttmp
        count += 1


U = np.concatenate((U1, U2, U3), axis=0)
V = np.concatenate((V1, V2, V3), axis=0)
T = np.concatenate((T1, T2, T3), axis=0)
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

