#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  NBUC_FM.py -- Building NBUC Feature Model
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Jan, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from roms_setup import near 
import scipy.io as sp
from mpl_toolkits.basemap import Basemap 
import seawater.csiro as sw
import netCDF4 as nc
import sympy.solvers as sym
from sympy.abc import x, y

# FUNCTIONS

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

# Basic Settings ##################################################################

lvel     = (-0.2, 0.2)
dataplot = 'n'

# LOOKING AT DATA
NbucSouth = sp.loadmat('NbucSouth.mat')
lonS = NbucSouth['lonsouth'][:]
zS   = NbucSouth['zsouth'][:]
vS   = NbucSouth['vsouth'][:]
vS   = np.ma.masked_where(vS > 10, vS)
vST  = np.ma.masked_where(vS < 0, vS)
vST  = np.ma.masked_where(lonS < -38, vST)
TS   = transp(lonS, zS, vST); strS = str( TS.round(2) ) +' Sv'

NbucCenter = sp.loadmat('NbucCenter.mat')
lonC = NbucCenter['loncenter'][:]
zC   = NbucCenter['zcenter'][:]
vC   = NbucCenter['vcenter'][:]
vC   = np.ma.masked_where(vC > 10, vC)
vCT  = np.ma.masked_where(vC < 0, vC)
vCT  = np.ma.masked_where(lonC < -38, vCT)
TC   = transp(lonC, zC, vCT); strC = str( TC.round(2) ) +' Sv'

NbucNorth = sp.loadmat('NbucNorth.mat')
lonN = NbucNorth['lonnorth'][:]
zN   = NbucNorth['znorth'][:]
vN   = NbucNorth['vnorth'][:]
vN   = np.ma.masked_where(vN > 10, vN)
vNT  = np.ma.masked_where(vN < 0, vN)
vNT  = np.ma.masked_where(lonN < -37.5, vNT)
TN   = transp(lonN, zN, vNT); strN = str( TN.round(2) ) +' Sv'

if dataplot == 'y':

	plt.figure(1,figsize=(12,12),facecolor='w')

	p1 = plt.subplot(3, 1, 1)
	p1.contourf(lonN, zN, vN, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
	plt.axis([lonN.min(), lonN.max(), -1000, 0])
	plt.colorbar()
	plt.text(lonN.mean(), -900, strN)
	p1.set_title("NBUC North")
	#p1.set_xlabel("Latitude")
	p1.set_ylabel("Z [m]")

	p2 = plt.subplot(3, 1, 2)
	p2.contourf(lonC, zC, vC, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
	plt.colorbar()
	plt.text(lonC.mean(), -900, strC)
	plt.axis([lonC.min(), lonC.max(), -1000, 0])
	p2.set_title("NBUC Center")
	# p2.set_xlabel("Latitude")
	p2.set_ylabel("Z [m]")

	p3 = plt.subplot(3, 1, 3)
	p3.contourf(lonS, zS, vS, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
	plt.colorbar()
	plt.text(lonS.mean(), -900, strS)
	plt.axis([lonS.min(), lonS.max(), -1000, 0])
	p3.set_title("NBUC South")
	p3.set_xlabel("Longitude")
	p3.set_ylabel("Z [m]")

	plt.show()


# ======================================================
# CREATING ISOBATH-FOLLOWING NBUC FEATURE MODEL:
# ======================================================

# ======================================================
# loading roms grid to get limits and topography
print ' '
print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
print ' '
# I need a bigger grid to get the isobath
grdfile  = nc.Dataset('/home/rsoutelino/rsoutelino/myroms/phd_run/phd1_grd.nc')

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h     = grdfile.variables['h'][:]

# getting an isobath
plt.figure(); con = plt.contour(lonr, latr, h, levels=[100] )
col = con.collections[0]; paths = col.get_paths()
path0 = paths[0]; isob = path0.vertices; plt.close('all')

# now I load original small grid
grdfile  = nc.Dataset('/home/rsoutelino/rsoutelino/myroms/phd_run/phd8_grd.nc')

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h     = grdfile.variables['h'][:]

# ======================================================
# defining some parameters to create feature model
D     = 2 # length of the transect (or total length of the modeled jet) [degrees]
tgmax = 0.2 # tan of the maximum angle that the transect is allowed to have to the horizontal
dr    = 0.01 / D # horizontal spacing for the feature model (adimensional) [degree / length]
r  = D * np.arange(0, 1, dr) # adimensional horizontal transects
r0 = (1.0/3.0) * D # defining transect center

# ======================================================
# defining domain, buffering variables
looprange = range(0, len(isob), 20) # scanning the isobath
li = np.size(looprange); lj = np.size(r)
X  = np.zeros( [li , lj]); Y = X.copy()
dz = 1.0 # DO NOT CHANGE THAT!!!! If changed, require changes in the code bellow here
Z  = np.arange(-1500.0, 0.0+dz, dz); lk = np.size(Z)

U = np.zeros( [ lk , li , lj ] )
V = U.copy(); VS = U.copy()
v = np.zeros( [ lk , lj ] )

# ======================================================
# defining jet cross-slope structure v = v(x)

# jet width: let's try constant
d  = ( (150.0/2.0) / 111.0 ) / D # km converted to degrees, normalized by transect length 

# ======================================================
# defining velocity-axis vertical structure v0 = v0(y,z)

# Y-dependance:  
# Jet core depth will be a linear function from south to north, coming from 350 m (South)
# to 200 m (North)
# Core velocity will also be a linear function from south to north, to allow the transport
# to increase (this has to be checked with data!)
z0max = np.linspace(-600, -200, li) # keeping constant now for simplicity (Jan-2011)
v0max = np.linspace(0.5, 0.5, li) # keeping constant now for simplicity (Jan-2011)
v0 = np.zeros( [lk, li] )

# Z-dependance:
# NBUC core in Z0 m, with 0.14 m/s, decaying in a gauss curve until 0 m/s at bottom
#    this will also by Y-dependant, to allow increasing of the jet thickness
# another gaussian will be adopted to decay velocities to surface 
d1      = 100 / dz # gaussian vertical upper width, normalized by dz
d2      = np.linspace(360/dz, 360/dz, li) # gaussian lower width, normalized by dz

# starting the looping to create the slope-following NBUC-FM
print ' '
print '======== CREATING SLOPE-FOLLOWING NBUC-FM ========'
print ' '
i = -1 # initializind index counter
for c in looprange:
    print '    Transect ' + str(i) + ' / ' + str(np.size(looprange))
    i = i + 1
    x0 = isob[c:c+6, 0]; y0 = isob[c:c+6, 1]
    tgr, b = np.polyfit(x0, y0, 1) # finding isobath-tangent straight line
    x0 = x0[0]; y0 = y0[0]
    tgr = -1.0 / tgr; b = y0 - tgr * x0 # finding normal straight line
    if tgr >= tgmax: # enforcing maximun angle
       tgr = tgmax
    elif tgr <= -tgmax:
       tgr = -tgmax
    
    # =======================================
    # assembling vertical jet-core structure
    
    # upper part
    z1 = Z[z0max[i]/dz:]
    v01 = v0max[i] * np.exp( (-1)* (z1 - (z0max[i]/dz))**2 / (2*d1**2) )
    v0[z0max[i]/dz:, i] = v01

    # lower part
    z2  = Z[:z0max[i]/dz]
    v02 = v0max[i] * np.exp( (-1)* (z2 - (z0max[i]/dz))**2 / (2*d2[i]**2) )
    v0[:z0max[i]/dz, i] = v02
    
    # ==========================================================
    # writing NBUC-FM to transect based on cross-slope structure
    
    for k in range( 0, lk):
        v[k, :] = v0[k, i] * np.exp( (-1)* ( ( r-r0 )**2 / ( 2*d**2 ) ) )

    # georeferencing velocities and coordinates
    angr = np.arctan(tgr)   
    cosr, sinr  = np.cos(angr), np.sin(angr)
    X[i, :]     = r * cosr + x0
    Y[i, :]     = r * sinr + y0
    U[:, i, :]  = v * sinr * (-1)
    V[:, i, :]  = v * cosr
    VS[:, i, :] = v
    
#  ==========================================================
# COMPUTING GEOSTROPHICALLY BALANCED STATE VARIABLES
print ' '
print '======== COMPUTING GEOSTROPHICALLY BALANCED STATE VARIABLES ========'
print ' '
# integration the thermal wind equation:
# rho(x,z) = rho0(z) - rho_bar.f/g * int_0^L{ dv/dz dx}

# obtaining rho0 and rho_bar from WOA2009:
MeanDens = sp.loadmat('MeanDens.mat')
rho0  = MeanDens['potdens'][:].ravel(); rho0 = rho0[::-1]
zrho  = MeanDens['z'][:].ravel(); zrho = zrho[::-1]
salt0 = MeanDens['salt'][:].ravel(); salt0 = salt0[::-1] 
rho0  = np.interp(Z, zrho, rho0)
salt0 = np.interp(Z, zrho, salt0)
rho0.shape = (np.size(rho0), 1)
salt0.shape = (np.size(salt0), 1)
rho_bar = rho0.mean()

print '    Density'
# obtaining dv/dz:
dvdz = np.zeros([ lk , li , lj])
for i in range(0, li):
    vaux = np.squeeze(VS[:,i,:]) 
    aux  = np.array(np.gradient(vaux))
    dvdz[:,i,:] = np.squeeze(aux[0, :, :])

# obtaining dS [m]: where S is the cross-slope axis
S  = r * 111 * 1000; S, tmp = np.meshgrid(S, Z)
dS = np.array(np.gradient(S))
dS = np.squeeze(dS[1, :, :])

# constants
g = 9.8
f0 = 4 * np.pi * np.sin( np.pi * latr.mean()/180 ) / ( 24*3600 ) 

# COMPUTING DENSITY:
RHO = np.zeros([ lk , li , lj])
for i in range(0, li):
    aux    =  dvdz[:,i,:]
    rhoaux  = rho0 - ((rho_bar*f0) / g) * np.cumsum( aux*dS, axis=1 )
    RHO[:,i,:]  = rhoaux

# COMPUTING TEMPERATURE AND SALINITY
# linearized equation of seawater state
alpha = 2.2e-4
beta  = 8e-4
S0    = 37
rho_ref  = 1000.7

TEMP = np.zeros([ lk , li , lj])
SALT = np.zeros([ lk , li , lj])
print '    Temperature, Salinity'
for i in range(0, li):
    TEMP[:,i,:]  = ( -RHO[:,i,:]/rho_ref + 1 + beta*salt0 ) / alpha
    SALT[:,i,:]  = salt0 + 0.01 * TEMP[:,i,:]
	
# =============================================================================
# DEFINING ENCOMPASSING DOMAIN BASED ON ROMS GRID TO INTERPOLATE FEATURE MODEL:
# =============================================================================

del v
dx, dy  = 0.2, 0.2 # horizontal resolution of the new grid [degrees]
newdz   = 10        # vertical resolution of the new grid [m]

# creating lon, lat arrays based on ROMS grid plus extra 1 degree on each boundary
lon = np.arange( lonr.min()-1, lonr.max()+1, dx )
lat = np.arange( latr.min()-1, latr.max()+1, dy )

# create metric arrays for further transport computing
xsize = ( (lonr.max()+1) - (lonr.min()-1) ) * 111.0
ysize = ( (latr.max()+1) - (latr.min()-1) ) * 111.0

x = np.linspace( 0.0, xsize, np.size(lon) )
y = np.linspace( 0.0, ysize, np.size(lat) )
z   = np.arange(-1500.0, 0.0+newdz, newdz)

lonsec, tmp  = np.meshgrid(lon, z)
lon, lat = np.meshgrid(lon, lat)
xx, zz = np.meshgrid(x*1000, z)

# buffering velocity arrays:
li, lj = lon.shape; lk = np.size(z)
u = np.zeros( [ lk , li , lj ] ); v = u.copy(); 
temp = u.copy(); salt = u.copy();

# interpolating FM to new grid
print ' '
print '======== GRIDDING NBUC-FM ========'
print ' '
for k in range(0, lk-1):
    print '    U:     Z Level = ' + str(-1*z[k]) + ' m'
    u[k,...] = griddata( X.ravel(), Y.ravel(), U[z[k],...].ravel(), lon, lat )
    print '    V:             = '
    v[k,...] = griddata( X.ravel(), Y.ravel(), V[z[k],...].ravel(), lon, lat )
    print '    TEMP:          = '
    temp[k,...] = griddata( X.ravel(), Y.ravel(), TEMP[z[k],...].ravel(), lon, lat )
    print '    SALT:          = '
    salt[k,...] = griddata( X.ravel(), Y.ravel(), SALT[z[k],...].ravel(), lon, lat )
    print ' '; print ' ';

# removing  nans
u[np.isnan(u)] = 0; v[np.isnan(v)] = 0
for k in range(0,lk-1):
    aux = temp[k,...]
    aux[np.isnan(aux)] = TEMP[z[k],-1,-1]
    temp[k,...] = aux
    aux = salt[k,...]
    aux[np.isnan(aux)] = SALT[z[k],-1,-1]
    salt[k,...] = aux
    
# ==================================
# PLOTTING AND COMPARING WITH DATA
# ==================================

# loading data
SecRomsWoa = sp.loadmat('SecRomsWoa.mat')
latdata = SecRomsWoa['latRomsWoa'][:]
zdata   = SecRomsWoa['zRomsWoa'][:]
udata   = SecRomsWoa['uRomsWoa'][:]
udata   = np.ma.masked_where(udata > 10, udata)

# computing transports
FMtranspS = transp(lonsec, zz, np.squeeze(v[:,10])) # to fit inside ROMS domain [hardcoded]
FMtranspN = transp(lonsec, zz, np.squeeze(v[:,-10])) # to fit inside ROMS domain [hardcoded]
str2 = str( FMtranspS.round(2) ) +' Sv'
str3 = str( FMtranspN.round(2) ) +' Sv'

# plotting
fig1 = plt.figure(1,figsize=(15,8),facecolor='w')

p1 = plt.subplot(2, 2, 1)
c1 = p1.contourf(lonS, zS, vST, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.axis([lonr.min(), lonr.max(), -1500, 0])
# plt.colorbar()
p1.set_xticklabels([])
plt.text(lon.mean(), -900, strS)
p1.set_title("NBUC OEII ROMS run: South")
# p1.set_xlabel("Latitude")
p1.set_ylabel("Z [m]")

p2 = plt.subplot(2, 2, 2)
p2.contourf(lonN, zN, vNT, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.axis([lonr.min(), lonr.max(), -1500, 0])
# plt.colorbar()
p2.set_xticklabels([])
p2.set_yticklabels([])
plt.text(lon.mean(), -900, strN)
p2.set_title("NBUC OEII ROMS run: North")
# p2.set_xlabel("Latitude")
# p2.set_ylabel("Z [m]")

p3 = plt.subplot(2, 2, 3)
p3.contourf(lonsec, zz, np.squeeze(v[:,10]), np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
# plt.colorbar()
plt.text(lon.mean(), -900, str2)
plt.axis([lonr.min(), lonr.max(), -1500, 0])
p3.set_title("NBUC FM: South")
p3.set_xlabel("Longitude")
p3.set_ylabel("Z [m]")

p4 = plt.subplot(2, 2, 4)
p4.contourf(lonsec, zz, np.squeeze(v[:,-10]), np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
# plt.colorbar()
p4.set_yticklabels([])
plt.text(lon.mean(), -900, str3)
plt.axis([lonr.min(), lonr.max(), -1500, 0])
p4.set_title("NBUC FM: North")
p4.set_xlabel("Longiude")
# p4.set_ylabel("Z [m]")

ax = fig1.add_axes([0.91, 0.1, 0.01, 0.8])
cbar = plt.colorbar(c1, cax=ax, orientation='vertical')
cbar.set_label('[m s$^{-1}$]')

plt.show()


# Plotting Horizontal slices of the NBUC FM

zb    = np.ma.masked_where(h >= 100, h)

m = Basemap(projection='merc',
     llcrnrlat = latr.min(), urcrnrlat = latr.max(),
     llcrnrlon = lonr.min(), urcrnrlon = lonr.max(),
     lat_ts=0, resolution='l')

lonrp, latrp = m(lonr, latr)

ss   = 12  
sc   = 4

lonp = lonsec[0,:]
#lonp, latp = np.meshgrid(lonp, lat)
#lonp = lonp[:, ::ss]
#latp = latp[:, ::ss]
lonp, latp = m(lon, lat)

u1 = np.squeeze( u[-1, :,:] );  v1 = np.squeeze( v[-1, :,:] )
u2 = np.squeeze( u[-70, :,:] ); v2 = np.squeeze( v[-70, :,:] )
#v1 = v1[:, ::ss]
#v2 = v2[:, ::ss]

fig2 = plt.figure(2,figsize=(14,8),facecolor='w')

pp1 = plt.subplot(1, 2, 1)
m.contourf(lonrp, latrp, zb, colors=('0.7'), alpha=0.5)
m.quiver(lonp, latp, u1, v1, scale=sc)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max()+2, 2),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max()+2, 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
pp1.set_title('NBUC FM - Surface Field')


pp2 = plt.subplot(1, 2, 2)
m.contourf(lonrp, latrp, zb, colors=('0.7'), alpha=0.5)
m.quiver(lonp, latp, u2, v2, scale=sc)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max()+2, 2),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max()+2, 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
pp2.set_title('NBUC FM - 400 m Field')

plt.show()



# computing geostrophic velocity to double-check
#temp1 = squeeze(temp[:,17,:])
#salt1 = squeeze(salt[:,17,:])
#gp    = sw.gpan(salt1, temp1,-zz)
#dgpdx = np.array(np.gradient(gp))
#dgpdx = np.squeeze(dgpdx[1,:,:])
#dgpdx = dgpdx / 1000 # [km] --> [m]
#gvel  = -dgpdx / f0
#gvel  = gvel - gvel[-1,:]

#for i in range(0, y.size):
#	figure(figsize=(15,5),facecolor='w')
#	subplot(1,2,1)
#	contourf(lon, zz, np.squeeze(temp[:,i,:]), np.arange(7, 29, 1), cmap=plt.cm.Spectral_r)
#	colorbar()
#	
#	subplot(1,2,2)
#	contourf(lon, zz, np.squeeze(v[:,i,:]), np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
#	colorbar()
#	plt.savefig('vel-temp_'+ str(i+1)+'.png')
#
#for k in range(1, 10):
#	execstr =  'mv vel-temp_'+ str(k) +'.png frame_0'+ str(k) +'.png'
#	os.system(execstr)
#
#for k in range(10, 100):
#	execstr =  'mv vel-temp_'+ str(k) +'.png frame_'+ str(k) +'.png'
#	os.system(execstr)
#
#
#os.system('convert -loop 0 -delay 60 frame_??.png NBUC_fm_movie.gif')
#

# PREPARING FIELDS TO BUILD ROMS INITIAL FIELDS

# Subsampling feature model to avoid unecessary huge computations
#ix = 10
#iy = 1
#iz = 10

#long = lon[::ix]
#latg = lat[::iy]
#zg   =    z[::iz]
#temp = temp[::iz, ::iy, ::ix]
#salt = salt[::iz, ::iy, ::ix]
#v    =    v[::iz, ::iy, ::ix]

# Expanding to complete the domain i.e. (buffer zones)
#     repeat values in the boundaries (N, S, E, W, Z)


#longW = ([long[0]-3, long[0]-2, long[0]-1])    # adding 3 degree west
#longE = ([long[-1]+1, long[-1]+2, long[-1]+3]) # adding 3 degree east
#long  = np.hstack((longW, long, longE))        # updating long

#latgS = ([latg[0]-3, latg[0]-2, latg[0]-1])    # adding 3 degree south
#latgN = ([latg[-1]+1, latg[-1]+2, latg[-1]+3]) # adding 3 degree north
#latg  = np.hstack((latgS, latg, latgN))        # updating latg

#zgB   = ([ -5500, -5000 , -4000, -3000, -2000 ]) # adding coordinate until bottom
#zg    = np.hstack((zgB, zg))                   # updating zg

## expanding velocity
#a, b, c = v.shape
#vW      = v[:,:,0];  vW.shape = (a,b,1); vW = vW.repeat(3, axis=2)
#vE      = v[:,:,-1]; vE.shape = (a,b,1); vE = vE.repeat(3, axis=2)
#v       = np.concatenate((vW, v, vE), axis=2); a, b, c = v.shape
#vS      = v[:,0,:];  vS.shape = (a,1,c); vS = vS.repeat(3, axis=1)
#vN      = v[:,-1,:]; vN.shape = (a,1,c); vN = vN.repeat(3, axis=1)
#v       = np.concatenate((vS, v, vN), axis=1); a, b, c = v.shape
#vB      = v[0,:,:];  vB.shape = (1,b,c); vB = vB.repeat(5, axis=0)
#v       = np.concatenate((vB, v), axis=0)

## expanding salinity
#a, b, c = salt.shape
#saltW   = salt[:,:,0];  saltW.shape = (a,b,1); saltW = saltW.repeat(3, axis=2)
#saltE   = salt[:,:,-1]; saltE.shape = (a,b,1); saltE = saltE.repeat(3, axis=2)
#salt    = np.concatenate((saltW, salt, saltE), axis=2); a, b, c = salt.shape
#saltS   = salt[:,0,:];  saltS.shape = (a,1,c); saltS = saltS.repeat(3, axis=1)
#saltN   = salt[:,-1,:]; saltN.shape = (a,1,c); saltN = saltN.repeat(3, axis=1)
#salt    = np.concatenate((saltS, salt, saltN), axis=1);  a, b, c = salt.shape
#saltB   = salt[0,:,:];  saltB.shape = (1,b,c); saltB = saltB.repeat(5, axis=0)
#salt    = np.concatenate((saltB, salt), axis=0)

## expanding temperature
#a, b, c = temp.shape
#tempW   = temp[:,:,0];  tempW.shape = (a,b,1); tempW = tempW.repeat(3, axis=2)
#tempE   = temp[:,:,-1]; tempE.shape = (a,b,1); tempE = tempE.repeat(3, axis=2)
#temp    = np.concatenate((tempW, temp, tempE), axis=2); a, b, c = temp.shape
#tempS   = temp[:,0,:];  tempS.shape = (a,1,c); tempS = tempS.repeat(3, axis=1)
#tempN   = temp[:,-1,:]; tempN.shape = (a,1,c); tempN = tempN.repeat(3, axis=1)
#temp    = np.concatenate((tempS, temp, tempN), axis=1);  a, b, c = temp.shape
#tempB   = temp[0,:,:];  tempB.shape = (1,b,c); tempB = tempB.repeat(5, axis=0)
#temp    = np.concatenate((tempB, temp), axis=0)


matdict = {'lon':lon, 'lat': lat, 'z': z, 'temp':temp, 'salt':salt, 'u':u, 'v':v}
sp.savemat('NBUC_FM.mat', matdict)

# creating and saving barotropic NBUC-FM
ubar = u.mean(axis=0)
vbar = v.mean(axis=0)

matdict = {'lon':lon, 'lat': lat, 'ubar':ubar, 'vbar':vbar}
sp.savemat('NBUC-BT_FM.mat', matdict)

# creating and saving barotropic 10x stronger NBUC-FM
ubar = ubar * 10
vbar = vbar * 10

matdict = {'lon':lon, 'lat': lat, 'ubar':ubar, 'vbar':vbar}
sp.savemat('NBUC10x-BT_FM.mat', matdict)






