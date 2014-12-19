#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  sec_fm.py -- Building SEC Feature Model
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Sep, 2010
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from roms_setup import near 
import scipy.io as sp
import seawater.csiro as sw

# FUNCTIONS

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

	transp = np.abs(dy) * np.abs(dz) * u; transp = transp.sum()
	transp = transp / 1e6

	return transp

# Basic Settings ##################################################################

lvel    = (-0.10, 0.10)
latlim  = (-25.0, -5.0)

# DEFINING SECTIONAL DOMAIN:

secsize = ( latlim[1] - latlim[0] ) * 111.0

y   = np.arange(0.0, secsize, 5.0) # [km] 
z   = np.arange(-1500.0, 1.0, 1.0)
lat = np.linspace( latlim[0], latlim[1], np.size(y) )
yy, zz = np.meshgrid(y*1000, z) 

# DEFINING PARAMETERS:

# defining velocity-axis vertical structure v0 = v0(z)
# SEC core in 150 m, with 0.092 m/s, decaying in a gauss curve until 0 m/s at 900 m
# from surface to 150 m, another gaussian will be adopted to bring from 0.09 to 0.092 cm/s 

u0 = np.zeros( np.size(z) )

u0surf = 0.088
u0bot  = 0.0
u0max  = 0.092
z0max  = -150.0

d1      = 42
z1 = z[z0max:]
u01 = (u0max - u0surf) * np.exp( (-1)* (z1 - (z0max))**2 / (2*d1**2) ) + u0surf
u0[z0max:] = -u01

d2  = 225
z2  = z[:z0max]
u02 = (u0max - u0bot) * np.exp( (-1)* (z2 - (z0max))**2 / (2*d2**2) ) + u0bot
u0[:z0max] = -u02


# depth-dependent position of the jet axis: y0 = y0(z)
y0 = np.zeros( np.size(z) )

# SEC axis position remais the same from surface to core-depth
f2 = np.array(near(lat, -14.5)).min()
y0[z0max:] = y[f2]


# SEC axis position migrates gently to the south, from 14.5 S to 16 S
y2 = y0[:z0max]
f1 = np.array(near(lat, -16)).min()
f2 = np.array(near(lat, -14.5)).min()
y2 = np.linspace( y[f1], y[f2], np.size(y2) ) # linear decay with depth
y0[:z0max] = y2

# jet width: let's try constant
d  = 170
 
# CREATING THE JET

u = np.zeros([ np.size(z) , np.size(y) ])

for k in np.arange( 0, np.size(z) ):
	u[k] = u0[k] * np.exp( (-1)* ( ( y-y0[k] )**2 / ( 2*d**2 ) ) )


# PLOTTING AND COMPARING WITH DATA

# loading data
SecRomsWoa = sp.loadmat('SecRomsWoa.mat')
latdata = SecRomsWoa['latRomsWoa'][:]
zdata   = SecRomsWoa['zRomsWoa'][:]
udata   = SecRomsWoa['uRomsWoa'][:]
udata   = np.ma.masked_where(udata > 10, udata)

# computing transports
latsec, tmp  = np.meshgrid(lat, z)
udataT   = np.ma.masked_where(udata > 0, udata)
DATAtransp = transp(latdata, zdata, udataT)
str1 = str( DATAtransp.round(2) ) + ' Sv'
uT   = np.ma.masked_where(u > 0, u)
FMtransp = transp(latsec, zz, uT)
str2 = str( FMtransp.round(2) ) +' Sv'

# plotting
plt.figure(1,figsize=(12,9),facecolor='w')

p1 = plt.subplot(2, 1, 1)
plt.contourf(latdata, zdata, udata, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.axis([lat.min(), lat.max(), -1000, 0])
plt.colorbar()
plt.text(lat.mean(), -900, str1)
p1.set_title("SEC section - WOA2009 ROMS run")
#p1.set_xlabel("Latitude")
p1.set_ylabel("Z [m]")

p2 = plt.subplot(2, 1, 2)
plt.contourf(latsec, zz, u, np.arange( lvel[0], lvel[1]+0.01, 0.01), cmap=plt.cm.RdBu)
plt.colorbar()
plt.text(lat.mean(), -900, str2)
plt.axis([lat.min(), lat.max(), -1000, 0])
p2.set_title("SEC Feature Model")
p2.set_xlabel("Latitude")
p2.set_ylabel("Z [m]")

plt.show()

# COMPUTING GEOSTROPHICALLY BALANCED STATE VARIABLES

# integration the thermal wind equation:
# rho(y,z) = rho0(z) + rho_bar.f/g * int_0^L{ du/dz dy}

# obtaining rho0 and rho_bar from WOA2009:
MeanDens = sp.loadmat('MeanDens.mat')
rho0  = MeanDens['potdens'][:].ravel(); rho0 = rho0[::-1]
zrho  = MeanDens['z'][:].ravel(); zrho = zrho[::-1]
salt0 = MeanDens['salt'][:].ravel(); salt0 = salt0[::-1] 
rho0  = np.interp(z, zrho, rho0)
salt0 = np.interp(z, zrho, salt0)
rho0.shape = (np.size(rho0), 1)
salt0.shape = (np.size(salt0), 1)
rho_bar = rho0.mean()

# obtaining dv/dz:
dudz  = np.array(np.gradient(u))
dudz  = np.squeeze(dudz[0, :, :])

# obtaining dy:
dy = np.array(np.gradient(yy))
dy = np.squeeze(dy[1, :, :])

# constants
g = 9.8
f0 = 4 * np.pi * np.sin( np.pi * lat.mean()/180 ) / ( 24*3600 ) 

# COMPUTING DENSITY:
rho  = rho0 + ((rho_bar*f0) / g) * np.cumsum( dudz*dy, axis=1 )

# COMPUTING TEMPERATURE AND SALINITY
# linearized equation of seawater state
alpha = 2.2e-4
beta  = 8e-4
S0    = 37
rho_ref  = 1000.7

temp  = ( -rho/rho_ref + 1 + beta*salt0 ) / alpha
salt  = salt0 + 0.01 * temp

# computing geostrophic velocity to double-check
#gp    = sw.gpan(salt, temp,-zz)
#dgpdy = np.array(np.gradient(gp))
#dgpdy = np.squeeze(dgpdy[1,:,:])
#dgpdy = dgpdy / 1000 # [km] --> [m]
#gvel  = dgpdy / f0
#gvel  = gvel - gvel[0,:]


# PREPARING FIELDS TO BUILD ROMS INITIAL FIELDS

# Subsampling feature model to avoid unecessary huge computations
iy = 10
iz = 1

latg =  lat[::iy]
zg   =    z[::iz]
temp = temp[::iz, ::iy]
salt = salt[::iz, ::iy]
u    =    u[::iz, ::iy]

# Expanding to complete the domain i.e. (buffer zones)
#     repeat values in the boundaries (N, S, Z)

latgS = ([latg[0]-3, latg[0]-2, latg[0]-1])    # adding 3 degree south
latgN = ([latg[-1]+1, latg[-1]+2, latg[-1]+3]) # adding 3 degree north
latg  = np.hstack((latgS, latg, latgN))        # updating latg
zgB   = ([ -5500, -5000 , -4000, -3000, -2000 ]) # adding coordinate until bottom
zg    = np.hstack((zgB, zg))                   # updating zg

# expanding velocity
a, b = u.shape
uS   = u[:,0];  uS.shape = (a,1); uS = uS.repeat(3, axis=1)
uN   = u[:,-1]; uN.shape = (a,1); uN = uN.repeat(3, axis=1)
u    = np.concatenate((uS, u, uN), axis=1); a, b = u.shape
uB   = u[0,:];  uB.shape = (1,b); uB = uB.repeat(5, axis=0)
u    = np.concatenate((uB, u), axis=0)

# expanding salinity
a, b = salt.shape
saltS   = salt[:,0];  saltS.shape = (a,1); saltS = saltS.repeat(3, axis=1)
saltN   = salt[:,-1]; saltN.shape = (a,1); saltN = saltN.repeat(3, axis=1)
salt    = np.concatenate((saltS, salt, saltN), axis=1); a, b = salt.shape
saltB   = salt[0,:];  saltB.shape = (1,b); saltB = saltB.repeat(5, axis=0)
salt    = np.concatenate((saltB, salt), axis=0)

# expanding temperature
a, b = temp.shape
tempS   = temp[:,0];  tempS.shape = (a,1); tempS = tempS.repeat(3, axis=1)
tempN   = temp[:,-1]; tempN.shape = (a,1); tempN = tempN.repeat(3, axis=1)
temp    = np.concatenate((tempS, temp, tempN), axis=1); a, b = temp.shape
tempB   = temp[0,:];  tempB.shape = (1,b); tempB = tempB.repeat(5, axis=0)
temp    = np.concatenate((tempB, temp), axis=0)

matdict = {'lat': latg, 'z': zg, 'temp':temp, 'salt':salt, 'u':u}
sp.savemat('SEC_FM.mat', matdict)






