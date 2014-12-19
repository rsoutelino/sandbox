#!/usr/bin/env python
###################################################################################
#  Script  plot_SEC.py -- Visualization of SEC sectional flow fields
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
	i1   = int( near( aux, lonlim[0] )[0] ) 
	i2   = int( near( aux, lonlim[1] )[0] )
	aux  = lat[:,0] 
	j1   = int( near( aux, latlim[0] )[0] )
	j2   = int( near( aux, latlim[1] )[0] )

	return i1, i2, j1, j2

def transp(lat, z, u):
	"""
	Slices some 3D field within some lon, lat limits
	Be carefull with the aproximation on distance computing
	"""
	dy = np.diff(lat, axis=1) * 111 * 1000  # valid only for low latitudes!!!
	aux  = dy[:,0]; aux.shape = (np.size(aux), 1)
	dy = np.concatenate( (dy, aux), axis=1) 

	dz = np.diff(z, axis=0)
	aux  = dz[0,:]; aux.shape = (1, np.size(aux))
	dz = np.concatenate( (dz, aux), axis=0)

	transp = dy * dz * u; transp = transp.sum()
	transp = transp / 1e6

	return transp

# Basic settings ##################################################################

lvel    = (-0.10, 0.10)
zmax    = -1000         # treshold depth to compute volume transport
lonlim  = (-36.3, -34)    # SEC's flow is going to be an average of these meridians
latlim  = (-20, -12)

# OE2 data - from ROMS diag RUN ###################################################
# I'm doing time-average too (4 days)

oe2file = '/Users/rsoutelino/rsoutelino/myroms/leste2_run/leste2prog_avg.nc'
stabday = 4   # day in which mean KE stabilized
ncfile = nc.Dataset(oe2file)
grd   = nc.Dataset(oe2file[:-6] + 'grd.nc')
lon   = ncfile.variables['lon_u'][:]
lat   = ncfile.variables['lat_u'][:]
u     = ncfile.variables['u'][stabday:-1,...]; u = u.mean(axis=0)
u     = np.ma.masked_where(u > 10, u)
z     = get_depths(ncfile, grd, stabday-1, 'u')
k, i, j = u.shape

# slicing sections and computing transports

i1, i2, j1, j2 = slice(lon, lat, lonlim, latlim)
u   = np.squeeze(   u[:, j1:j2, i1:i2] ); u = u.mean(axis=2)
lat = np.squeeze(    lat[j1:j2, i1:i2] ); lat = lat.mean(axis=1)
z   = np.squeeze(   z[:, j1:j2, i1:i2] ); z = z.mean(axis=2)
lat.shape = ( 1, np.size(lat) )
lat = lat.repeat(k, axis=0)
u  = np.ma.masked_where(z < zmax, u)

Tv = transp(lat, z, u)

if Tv > 0:
	strTv = str( Tv.round(3) ) +' Sv: Outcoming'
else:
	strTv = str( Tv.round(3) ) +' Sv: Incoming'


plt.figure(1,figsize=(12,5),facecolor='w')
plt.contourf(lat, z, u, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.title('SEC section - OEII ROMS run')
plt.xlabel('Latitude')
plt.ylabel('Depth')
plt.text(lat.mean(), zmax + 100, strTv)
plt.colorbar()
plt.axis([lat.min(), lat.max(), zmax, 0])
plt.show()

# WOA2009 summer - from ROMS run ###################################################
# I'm doing time-average too (10 days)

woafile = '/Users/rsoutelino/rsoutelino/myroms/phd_run/phd2_avg.nc'
stabday = 10   # day in which mean KE stabilized
ncfile = nc.Dataset(woafile)
grd   = nc.Dataset(woafile[:-6] + 'grd.nc')
lon   = ncfile.variables['lon_u'][:]
lat   = ncfile.variables['lat_u'][:]
u     = ncfile.variables['u'][stabday:stabday+10,...]; u = u.mean(axis=0)
u     = np.ma.masked_where(u > 10, u)
z     = get_depths(ncfile, grd, stabday-1, 'u')
k, i, j = u.shape

# slicing sections and computing transports

i1, i2, j1, j2 = slice(lon, lat, lonlim, latlim)
u   = np.squeeze(   u[:, j1:j2, i1:i2] ); u = u.mean(axis=2)
lat = np.squeeze(    lat[j1:j2, i1:i2] ); lat = lat.mean(axis=1)
z   = np.squeeze(   z[:, j1:j2, i1:i2] ); z = z.mean(axis=2)
lat.shape = ( 1, np.size(lat) )
lat = lat.repeat(k, axis=0)
u  = np.ma.masked_where(z < zmax, u)

Tv = transp(lat, z, u)

if Tv > 0:
	strTv = str( Tv.round(3) ) +' Sv: Outcoming'
else:
	strTv = str( Tv.round(3) ) +' Sv: Incoming'


plt.figure(2,figsize=(12,5),facecolor='w')
plt.contourf(lat, z, u, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.title('SEC section - WOA2009 ROMS run')
plt.xlabel('Latitude')
plt.ylabel('Depth')
plt.text(lat.mean(), zmax + 100, strTv)
plt.colorbar()
plt.axis([lat.min(), lat.max(), zmax, 0])
plt.show()
