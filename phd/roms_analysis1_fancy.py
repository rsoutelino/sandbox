#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  roms_analysis1.py -- Analysis of ROMS fields
#  - 3-boundary transport time-series
#
#  Last Modification: May, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
from cookb_signalsmooth import smooth

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near

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

	transp = np.abs(dx) * np.abs(dz) * v; transp = transp.sum()
	transp = transp / 1e6

	return transp
	

# PYTHON SCRIPT START ######################################################

# ROMS SETTINGS ==================================================================

expt        = 'phd15'      # case nickname
filetype    = 'avg'        # could be his, avg, rst, ini, clm, ...
days        = range(0,90,1) # length of the series
vsc         = (-0.3, 0.3+0.03)
vst         = 0.03        # contour interval

niN         = -13         # latitude for northern section
niS         = -19         # latitude for southern section
njE         = -34         # longitude for eastern section

XlimN       = (-38, -35)  # x axis limits for northern section
XlimS       = (-39, -35)  # x axis limits for southern section
YlimE       = (-22, -10)  # x axis limits for eastern section

Zlim        = (-1500, 0)  # z axis limits

# READING ROMS FIELDS ============================================================

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
roms = run_setup(expt + '_run.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(roms.rundir + roms.run_name + '_grd.nc')

# assigning some variables from grid file
lon   = grdfile.variables['lon_rho'][:]
lat   = grdfile.variables['lat_rho'][:]
latu  = grdfile.variables['lat_u'][:]
latv  = grdfile.variables['lat_v'][:]
lonu  = grdfile.variables['lon_u'][:]
lonv  = grdfile.variables['lon_v'][:]
h     = grdfile.variables['h'][:]

print ' \n' + '==> ' + '  READING CHOSEN ROMS OUTPUT NETCDF FILE ...\n' + ' '
outfile  = netCDF4.Dataset(roms.rundir + roms.run_name + '_' + filetype + '.nc')

TSbc = []; TSnbuc = []; TN = []; TE = [];

for l in days:
    U    = outfile.variables['u'][l,...]
    V    = outfile.variables['v'][l,...]
    print " DAY = " + str(l+1)
    #km, im, jm = T.shape
    zu   = get_depths(outfile,grdfile,l,'u')
    zv   = get_depths(outfile,grdfile,l,'v')

    f    = latv[:,0]    
    fsec = near(f,niN)
    xN   = lonv[0,:]
    zN   = np.squeeze( zv[:, fsec, :] )
    vN   = np.squeeze( V[: ,fsec , :] ) 
    
    fsec = near(f,niS)
    xS   = lonv[0,:]
    zS   = np.squeeze( zv[:, fsec, :] )
    vS   = np.squeeze( V[: ,fsec , :] ) 
 
    f = lonu[0, :]
    fsec = near(f,njE)
    yE   = latu[:, 0]
    zE   = np.squeeze( zu[:, :, fsec] )
    uE   = np.squeeze( U[: ,: , fsec] ) 

    xS.shape = (1, xS.size)
    xN.shape = (1, xN.size)
    yE.shape = (1, yE.size)

    xS = np.repeat(xS, roms.klevels, axis=0)
    xN = np.repeat(xN, roms.klevels, axis=0)
    yE = np.repeat(yE, roms.klevels, axis=0)
    
    vS = np.ma.masked_where(vS > 1e10, vS)
    vN = np.ma.masked_where(vN > 1e10, vN)
    uE = np.ma.masked_where(uE > 1e10, uE)

    #fig1 = plt.figure(1, figsize=(12,6), facecolor='w')
    #ax = plt.axes(axisbg=('0.7'))
    #plt.contourf(xN,zN,vN, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu)
    #plt.colorbar(); plt.axis([XlimN[0], XlimN[1], Zlim[0], Zlim[1]])
    #plt.xlabel('Longitude')
    #plt.ylabel('Depth')
    #plt.title('Cross-section Velocity - Day: '+ str(l+1))

    
    #fig2 = plt.figure(2, figsize=(12,6), facecolor='w')
    #ax = plt.axes(axisbg=('0.7'))
    #plt.contourf(xS,zS,vS, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu)
    #plt.colorbar(); plt.axis([XlimS[0], XlimS[1], Zlim[0], Zlim[1]])
    #plt.xlabel('Longitude')
    #plt.ylabel('Depth')
    #plt.title('Cross-section Velocity - Day: '+ str(l+1))
    
    
    #fig3 = plt.figure(3, figsize=(12,6), facecolor='w')
    #ax = plt.axes(axisbg=('0.7'))
    #plt.contourf(yE,zE,uE, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu)
    #plt.colorbar(); plt.axis([YlimE[0], YlimE[1], Zlim[0], Zlim[1]])
    #plt.xlabel('Latitude')
    #plt.ylabel('Depth')
    #plt.title('Cross-section Velocity - Day: '+ str(l+1))
    
    # northern section transport
    vN = np.ma.masked_where(xN > XlimN[1], vN)
    vN = np.ma.masked_where(vN < 0, vN)
    tv = transp(xN, zN, vN)
    TN = np.hstack((TN, tv))
    
    # southern section transport
    BC   = np.ma.masked_where(xN > XlimS[1], vS)
    BC   = np.ma.masked_where(BC > 0, BC)
    tv   = transp(xS, zS, BC)
    TSbc = np.hstack((TSbc, tv))
    
    NBUC = np.ma.masked_where(xN > XlimS[1], vS)
    NBUC = np.ma.masked_where(NBUC < 0, NBUC)
    tv   = transp(xS, zS, NBUC)
    TSnbuc = np.hstack((TSnbuc, tv))
    
    # eastern section transport
    #uE = np.ma.masked_where(uE > 0, uE)
    tv = transp(yE, zE, uE)
    TE = np.hstack((TE, tv))
    
    
plt.figure(1, figsize=(12,6), facecolor='w')
p1 = plt.subplot(221)
plt.plot(days, TN, 'b', linewidth=1.5); plt.grid()
plt.title('Transport @ Northern Boundary', fontsize=11)
p1.set_xticklabels(' ')
p1.set_xlim([days[0], days[-1]])
p1.set_ylim([5, 40])
plt.ylabel('Transport [Sv]', fontsize=11)

p2 = plt.subplot(222)
plt.plot(days, TE, 'g', linewidth=1.5); plt.grid()
plt.title('Transport @ Eastern Boundary', fontsize=11)
p2.set_xticklabels(' ')
plt.ylabel('Transport [Sv]', fontsize=11)
#plt.xlabel('Days')
p2.set_xlim([days[0], days[-1]])
p2.set_ylim([-25, -5])

p3 = plt.subplot(223)
plt.plot(days, TSbc, 'r', linewidth=1.5, label='BC'); plt.grid()
plt.plot(days, TSnbuc, 'b', linewidth=1.5, label='NBUC'); plt.legend()
plt.title('Transport @ Southern Boundary', fontsize=11)
plt.xlabel('Days', fontsize=11)
plt.ylabel('Transport [Sv]', fontsize=11)
p3.set_xlim([days[0], days[-1]])
p3.set_ylim([-15, 25])

data = np.loadtxt('phd15_ek.out')

dt = 300

t = data[:,0]
t = (t*dt)  / 86400 
ke = data[:,5]
ke = smooth(ke, window_len=201, window='hanning')

p3 = plt.subplot(224)
plt.plot(t, ke, 'k', linewidth=1.5); plt.grid()
plt.title('Domain-averaged Kinetic Energy', fontsize=11)
plt.xlabel('Days', fontsize=11)
plt.ylabel('J m$^{-2}$', fontsize=11)
p3.set_xlim([days[0], days[-1]])
#p3.set_ylim([-15, 25])
    
plt.show()
plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/transp_inside2.pdf')
