#!/usr/bin/env python
###################################################################################
#  Script  plot_NBUC.py -- Visualization of NBUC sectional flow fields
#  As support for making the Feature Model
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Sep, 2010
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4 as nc

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near 

# FUNCTIONS ######################################################################

def slice(lon, lat, lonlim, latlim):
	"""
	Slices some 3D field within some lon, lat limits
	"""
	aux  = lon[0,:] 
	i1   = np.array( near( aux, lonlim[0] ) ).min() 
	i2   = np.array( near( aux, lonlim[1] ) ).min()
	aux  = lat[:,0] 
	j1   = np.array( near( aux, latlim[0] ) ).min()
	j2   = np.array( near( aux, latlim[1] ) ).min()

	return i1, i2, j1, j2

def transp(lon, z, v):
	"""
	Slices some 3D field within some lon, lat limits
	Be carefull with the aproximation on distance computing
	"""
	dx = np.diff(lon, axis=1) * 111 * 1000  # valid only for low latitudes!!!
	aux  = dx[:,0]; aux.shape = (np.size(aux), 1)
	dx = np.concatenate( (dx, aux), axis=1) 

	dz = np.diff(z, axis=0)
	aux  = dz[0,:]; aux.shape = (1, np.size(aux))
	dz = np.concatenate( (dz, aux), axis=0)

	transp = dx * dz * v; transp = transp.sum()
	transp = transp / 1e6

	return transp

# Basic settings ##################################################################

lvel    = (-0.25, 0.25)
zmax    = -1000         # treshold depth to compute volume transport
lonlim  = (-40, -36)    # 
latlim  = (-17, -14)    # NBUC's flow is going to be an average of these meridians

# OE2 data - from ROMS diag RUN ###################################################
# I'm doing time-average too (4 days)

oe2file = '/Users/rsoutelino/rsoutelino/myroms/leste2_run/leste2geo_avg.nc'
stabday = 4   # day in which mean KE stabilized
ncfile = nc.Dataset(oe2file)
grd   = nc.Dataset(oe2file[:-6] + 'grd.nc')
lon   = ncfile.variables['lon_v'][:]
lat   = ncfile.variables['lat_v'][:]
v     = ncfile.variables['v'][1:stabday,...];  v = v.mean(axis=0)
v     = np.ma.masked_where(v > 10, v)
z     = get_depths(ncfile, grd, stabday-1, 'v')
k, i, j = v.shape

# slicing sections and computing transports

i1, i2, j1, j2 = slice(lon, lat, lonlim, latlim)
v   = np.squeeze(   v[:, j1:j2, i1:i2] ); v = v.mean(axis=1)
lon = np.squeeze(    lon[j1:j2, i1:i2] ); lon = lon.mean(axis=0)
z   = np.squeeze(   z[:, j1:j2, i1:i2] ); z = z.mean(axis=1)
lon.shape = ( 1, np.size(lon) )
lon = lon.repeat(k, axis=0)
v  = np.ma.masked_where(z < zmax, v)

Tv = transp(lon, z, v)

if Tv > 0:
	strTv = str( Tv.round(2) ) +' Sv'
else:
	strTv = str( Tv.round(2) ) +' Sv'


plt.figure(1,figsize=(12,5),facecolor='w')
plt.contourf(lon, z, v, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.title('NBUC section - OEII ROMS run')
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.text(lon.mean(), zmax + 100, strTv)
plt.colorbar()
plt.axis([lon.min(), lon.max(), zmax, 0])
plt.show()

# WOA2009 summer - from ROMS run ###################################################
# I'm doing time-average too (10 days)
STOP
lonlim  = (-40, -30)
latlim  = (-19, -18)  

woafile = '/Users/rsoutelino/rsoutelino/myroms/phd_run/phd2_avg.nc'
stabday = 10   # day in which mean KE stabilized
ncfile = nc.Dataset(woafile)
grd   = nc.Dataset(woafile[:-6] + 'grd.nc')
lon   = ncfile.variables['lon_v'][:]
lat   = ncfile.variables['lat_v'][:]
v     = ncfile.variables['v'][0:stabday,...]; v = v.mean(axis=0)
v     = np.ma.masked_where(v > 10, v)
z     = get_depths(ncfile, grd, stabday, 'v')
k, i, j = v.shape

# slicing sections and computing transports

i1, i2, j1, j2 = slice(lon, lat, lonlim, latlim)
v   = np.squeeze(   v[:, j1:j2, i1:i2] ); v = v.mean(axis=1)
lon = np.squeeze(    lon[j1:j2, i1:i2] ); lon = lon.mean(axis=0)
z   = np.squeeze(   z[:, j1:j2, i1:i2] ); z = z.mean(axis=1)
lon.shape = ( 1, np.size(lon) )
lon = lon.repeat(k, axis=0)
v  = np.ma.masked_where(z < zmax, v)

Tv = transp(lon, z, v)

if Tv > 0:
	strTv = str( Tv.round(2) ) +' Sv'
else:
	strTv = str( Tv.round(2) ) +' Sv'


plt.figure(2,figsize=(12,5),facecolor='w')
plt.contourf(lon, z, v, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.title('NBUC section - WOA2009 ROMS run')
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.text(lon.mean(), zmax + 100, strTv)
plt.colorbar()
plt.axis([lon.min(), lon.max(), zmax, 0])
plt.show()
