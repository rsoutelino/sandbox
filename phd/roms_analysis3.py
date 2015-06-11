#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  roms_analysis1.py -- Analysis of ROMS fields
#  - vertical section mesh-up, averages, etc
#
#  Last Modification: May, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4

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
days        = range(0,359) # length of the series


tsc         = (3, 30)    # color scale limits for temperature
tst         = 1.5          # contour interval
ssc         = (34.4, 37)   # color scale limits for salinity
sst         = (0.2)     # contour interval
vsc         = (-0.4, 0.03+0.4, 0.03)
vst         = 0.03        # contour interval

niN         = -13         # latitude for NBUC section
niE         = -34         # latitude for RC Eddy
niS         = -19         # longitude for Abrolhos Eddy

XlimN       = (-39, -34)  # x axis limits for northern section
XlimE       = (-23, -10)  # x axis limits for southern section
XlimS       = (-39, -34)  # x axis limits for eastern section

Zlim        = (-1000, 0)  # z axis limits

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

print ' \n' + '==> ' + '  COMPUTING TIME-AVERAGE ...\n' + ' '
T=0; S=0; U=0; V=0; W=0;
for l in days:
    T    = T + outfile.variables['temp'][l,...]
    S    = S + outfile.variables['salt'][l,...]
    U    = U + outfile.variables['u'][l,...]
    V    = V + outfile.variables['v'][l,...]

T = T / len(days); S = S / len(days); 
U = U / len(days); V = V / len(days); 

km, im, jm = T.shape

print ' \n' + '==> ' + '  GETTING DEPTHS OF S-LEVELS ...\n' + ' '
zt   = get_depths(outfile,grdfile,days[0],'temp')
zu   = get_depths(outfile,grdfile,days[0],'u')
zv   = get_depths(outfile,grdfile,days[0],'v')


# vertical ==============================================

f    = latv[:,0]    
fsec = near(f,niN)
xN   = lonv[0,:]
zN   = np.squeeze( zv[:, fsec, :] )
vN   = np.squeeze( V[: ,fsec , :] ) 
tN   = np.squeeze( T[: ,fsec , :] )
sN   = np.squeeze( S[: ,fsec , :] )
ztN  = np.squeeze( zt[: ,fsec , :] )

fsec = near(f,niS)
xS   = lonv[0,:]
zS   = np.squeeze( zv[:, fsec, :] )
vS   = np.squeeze( V[: ,fsec , :] ) 
tS   = np.squeeze( T[: ,fsec , :] )
sS   = np.squeeze( S[: ,fsec , :] )
ztS  = np.squeeze( zt[: ,fsec , :] )

f    = lonu[0,:] 
fsec = near(niE,f)
xE   = latu[:,0]
zE   = np.squeeze( zu[:, :, fsec] )
uE   = np.squeeze( U[:,  :, fsec] ) 
tE   = np.squeeze( T[: , :, fsec] )
sE   = np.squeeze( S[: , :, fsec] )
ztE  = np.squeeze( zt[: ,:, fsec] )

xS.shape = (1, xS.size)
xN.shape = (1, xN.size)
xE.shape = (1, xE.size)

xS = np.repeat(xS, roms.klevels, axis=0)
xN = np.repeat(xN, roms.klevels, axis=0)
xE = np.repeat(xE, roms.klevels, axis=0)

vS = np.ma.masked_where(vS > 1e10, vS)
vN = np.ma.masked_where(vN > 1e10, vN)
uE = np.ma.masked_where(uE > 1e10, uE)

### NORTH ------------------------------------------------------------------------------

fig1 = plt.figure(1, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xN,ztN,tN, np.arange(tsc[0], tsc[1], tst),cmap=plt.cm.Spectral_r, extend='both')
plt.colorbar(); plt.axis([XlimN[0], XlimN[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Temperature  @ '+ str(np.abs(niN)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niN)) +'S.png')
plt.contour(xN,zN,vN,10, colors='k')
plt.savefig('figures/'+ expt +'/vel-temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niN)) +'S.png')
   
fig2 = plt.figure(2, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xN,ztN,sN, np.arange(ssc[0], ssc[1], sst),cmap=plt.cm.YlOrBr, extend='both')
plt.colorbar(); plt.axis([XlimN[0], XlimN[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Salinity  @ '+ str(np.abs(niN)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/salt_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niN)) +'S.png')

fig3 = plt.figure(3, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xN,zN,vN, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimN[0], XlimN[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niN)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/vel_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niN)) +'S.png')


### EAST -----------------------------------------------------------------------------------

fig4 = plt.figure(4, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xE,ztE,tE, np.arange(tsc[0], tsc[1], tst),cmap=plt.cm.Spectral_r, extend='both')
plt.colorbar(); plt.axis([XlimE[0], XlimE[1], Zlim[0], Zlim[1]])
plt.xlabel('Latitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Temperature  @ '+ str(np.abs(niE)) +'$^\circ$W')
plt.savefig('figures/'+ expt +'/temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niE)) +'W.png')
plt.contour(xE,zE,uE,10, colors='k')
plt.savefig('figures/'+ expt +'/vel-temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niE)) +'W.png')
   
fig5 = plt.figure(5, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xE,ztE,sE, np.arange(ssc[0], ssc[1], sst),cmap=plt.cm.YlOrBr, extend='both')
plt.colorbar(); plt.axis([XlimE[0], XlimE[1], Zlim[0], Zlim[1]])
plt.xlabel('Latitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Salinity  @ '+ str(np.abs(niE)) +'$^\circ$W')
plt.savefig('figures/'+ expt +'/salt_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niE)) +'W.png')   
   
fig6 = plt.figure(6, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xE,zE,uE, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimE[0], XlimE[1], Zlim[0], Zlim[1]])
plt.xlabel('Latitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niE)) +'$^\circ$W')
plt.savefig('figures/'+ expt +'/vel_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niE)) +'W.png')
   
   
### SOUTH -----------------------------------------------------------------------------------

fig7 = plt.figure(7, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xS,ztS,tS, np.arange(tsc[0], tsc[1], tst),cmap=plt.cm.Spectral_r, extend='both')
plt.colorbar(); plt.axis([XlimS[0], XlimS[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Temperature  @ '+ str(np.abs(niS)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niS)) +'S.png')
plt.contour(xS,zS,vS,10, colors='k')
plt.savefig('figures/'+ expt +'/vel-temp_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niS)) +'S.png')
   
fig8 = plt.figure(8, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xS,ztS,sS, np.arange(ssc[0], ssc[1], sst),cmap=plt.cm.YlOrBr, extend='both')
plt.colorbar(); plt.axis([XlimS[0], XlimS[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Salinity  @ '+ str(np.abs(niS)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/salt_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niS)) +'S.png')   
     
fig9 = plt.figure(9, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xS,zS,vS, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimS[0], XlimS[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niS)) +'$^\circ$S')
plt.savefig('figures/'+ expt +'/vel_sec_day'+
   str(days[0]+1) +'-'+ str(days[-1]+1) +'_'+ str(np.abs(niS)) +'S.png')

plt.show()
