#!/usr/bin/env python
###################################################################################
#  Script  volcons.py -- Checking volume conservation in Brazil East Box
#  
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Sep, 2010
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import find
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import seawater as sw
from roms_setup import run_setup, get_depths, zlevs, near 

def distance(lat1, lon1, lat2, lon2):
	"""
	Caclulate distance between two lat lons in NM
	"""
	rad = np.pi / 180.0
	yDistance = (lat2 - lat1) * 60.00721
	xDistance = (np.cos(lat1 * rad) + np.cos(lat2 * rad)) *\
                (lon2 - lon1) * (60.10793 / 2)
	distance = np.sqrt( yDistance**2 + xDistance**2 )
	return distance * 1.15078


# Basic settings ##################################################################

oe2file = '/Users/rsoutelino/rsoutelino/myroms/leste2_run/leste2prog_avg.nc'
stabday = 4   # day in which mean KE stabilized
lvel    = (-0.25, 0.25)
zmax    = -1000 # treshold depth to compute volume transport

# defining the working BOX:

lonlim = (-40, -34)
latlim = (-21, -11)
m = Basemap(projection='merc', llcrnrlat=latlim[0], urcrnrlat=latlim[1],
	llcrnrlon=lonlim[0], urcrnrlon=lonlim[1], lat_ts=0, resolution='l')


# SCRIPT START ###################################################################

oe2nc = nc.Dataset(oe2file)
grd   = nc.Dataset(oe2file[:-6] + 'grd.nc')
lon   = oe2nc.variables['lon_rho'][:]
lat   = oe2nc.variables['lat_rho'][:]
lonu  = oe2nc.variables['lon_u'][:]
latu  = oe2nc.variables['lat_u'][:]
lonv  = oe2nc.variables['lon_v'][:]
latv  = oe2nc.variables['lat_v'][:]

u   = oe2nc.variables['u'][stabday-1,]
v   = oe2nc.variables['v'][stabday-1,...]

# masking the arrays
u = np.ma.masked_where(u > 10, u)
v = np.ma.masked_where(v > 10, v)

ku, iu, ju = u.shape
kv, iv, jv = v.shape

zu   = get_depths(oe2nc, grd, stabday-1, 'u')
zv   = get_depths(oe2nc, grd, stabday-1, 'v')

# SLICING THE SECTIONS and COMPUTING TRANSPORTS

# South boundary
vS   = np.squeeze(   v[:, 2, :] )
lonS = np.squeeze(   lonv[2, :] )
zS   = np.squeeze(  zv[:, 2, :] )
lonS.shape = ( 1, np.size(lonS) )
lonS = lonS.repeat(ku, axis=0)
vSt  = np.ma.masked_where(zS < zmax, vS)

dx = np.diff(lonS, axis=1) * 111 * 1000
aux  = dx[:,0]; aux.shape = (np.size(aux), 1)
dx = np.concatenate( (dx, aux), axis=1) 

dz = np.diff(zS, axis=0)
aux  = dz[0,:]; aux.shape = (1, np.size(aux))
dz = np.concatenate( (dz, aux), axis=0)

transpS = dx * dz * vSt; transpS = transpS.sum()
transpS = transpS / 1e6


# East boundary
uE   = np.squeeze(   u[:, 2:-2, -3] )
latE = np.squeeze(   latu[2:-2, -3] )
zE   = np.squeeze(  zu[:, 2:-2, -3] )
latE.shape = ( 1, np.size(latE) )
latE = latE.repeat(ku, axis=0)
uEt  = np.ma.masked_where(zE < zmax, uE)

dy = np.diff(latE, axis=1) * 111 * 1000
aux  = dy[:,0]; aux.shape = (np.size(aux), 1)
dy = np.concatenate( (dy, aux), axis=1) 

dz = np.diff(zE, axis=0)
aux  = dz[0,:]; aux.shape = (1, np.size(aux))
dz = np.concatenate( (dz, aux), axis=0)

transpE = dy * dz * uEt; transpE = transpE.sum()
transpE = transpE / 1e6

# North boundary
vN   = np.squeeze(   v[:,-3, :] )
lonN = np.squeeze(  lonv[-3, :] )
zN   = np.squeeze(  zv[:,-3, :] )
lonN.shape = ( 1, np.size(lonN) )
lonN = lonN.repeat(ku, axis=0)
vNt  = np.ma.masked_where(zN < zmax, vN)

dx = np.diff(lonN, axis=1) * 111 * 1000
aux  = dx[:,0]; aux.shape = (np.size(aux), 1)
dx = np.concatenate( (dx, aux), axis=1) 

dz = np.diff(zN, axis=0)
aux  = dz[0,:]; aux.shape = (1, np.size(aux))
dz = np.concatenate( (dz, aux), axis=0)

transpN = dx * dz * vNt; transpN = transpN.sum()
transpN = transpN / 1e6

# Assigning transport direction: 
if transpS > 0:
	strS = str( transpS.round(2) ) +' Sv: Incoming'
else:
	strS = str( transpS.round(2) ) +' Sv: Outgoing'

if transpE < 0:
	strE = str( transpE.round(2) ) +' Sv: Incoming'
else:
	strE = str( transpE.round(2) ) +' Sv: Outgoing'

if transpN < 0:
	strN = str( transpN.round(2) ) +' Sv: Incoming'
else:
	strN = str( transpN.round(2) ) +' Sv: Outgoing'


# PLOTTING ##################################################################

fig1 = plt.figure(1,figsize=(15,9),facecolor='w')

p1 = plt.subplot(1, 3, 1)
p1.set_axis_bgcolor('0.5')
c1 = p1.pcolormesh(lonS, zS, vSt, vmin=lvel[0], vmax=lvel[1], cmap=plt.cm.RdBu)
p1.set_title('South boundary')
p1.set_xlabel('Longitude')
p1.set_ylabel('Depth')
p1.text(lonS.mean(), zS.mean(), strS)

p2 = plt.subplot(1, 3, 2)
p2.set_axis_bgcolor('0.5')
c2 = p2.pcolormesh(latE, zE, uEt, vmin=lvel[0], vmax=lvel[1], cmap=plt.cm.RdBu)
p2.set_title('East boundary')
p2.set_xlabel('Latitude')
p2.set_yticklabels(' ')
p2.text(latE.mean(), zE.mean(), strE)


p3 = plt.subplot(1, 3, 3)
p3.set_axis_bgcolor('0.5')
c3 = p3.pcolormesh(lonN, zN, vNt, vmin=lvel[0], vmax=lvel[1], cmap=plt.cm.RdBu)
p3.set_title('North boundary')
p3.set_xlabel('Longitude')
p3.set_yticklabels(' ')
p3.text(lonN.mean(), zN.mean(), strN)

ax = fig1.add_axes([0.91, 0.1, 0.01, 0.8])
cbar = plt.colorbar(c1, cax=ax, orientation='vertical')
cbar.set_label('[m s$^{-1}$]')


