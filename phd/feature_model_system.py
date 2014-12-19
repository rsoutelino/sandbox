#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#  Building 3-D Idealized velocity and Mass Fields for the BC-NBUC system
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Oct, 2012
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from roms_setup import near 
import scipy.io as sp
from mpl_toolkits.basemap import Basemap 
import seawater.csiro as sw
import netCDF4 as nc
from cookb_signalsmooth import smooth


# INPUT PARAMETERS #############################################################

# General Parameters ===========================================================

largegridfile = "/home/rsoutelino/myroms/phd_run/phd1_grd.nc"
smallgridfile = "/home/rsoutelino/myroms/phd_run/phd8_grd.nc"
climatology_profile = "MeanDens.mat"

# Common FM Parameters =========================================================

D         = 3.                    # length of the transect or total length of  
                                  # the modeled jet [degrees]
                                     
tgmax     = 0.15                  # tan of the maximum angle that a transect is
                                  #     allowed to have in respect to a parallel
                                     
dr        = 0.1                   # horizontal spacing [degree]
                                     
ln        = 30.                   # jump on the isobath to define the
                                  #     number of transects

zmax      = 1500.                 # maximum depth of the system

dz        = 1.                    # vertical resolution of the transects [m]
                                  #     DO NOT CHANGE THAT!!!! If changed, 
                                  #     require changes in the code

delta     = 100.                  # jets width [km]

lims      = (-60, -30, -24, -8)   # longitude / latidude limits

# NBUC-FM parameters ===========================================================

z0_NBUC_S = 500.                  # incoming NBUC core depth [m]

z0_NBUC_N = 200.                  # outgoing NBUC core depth [m]

v0_NBUC_S = 0.2                   # incoming NBUC core velocity [m/s]

v0_NBUC_N = 0.5                   # outgoing NBUC core velocity [m/s]

ds_NBUC   = 100.                  # NBUC thickness from core to surface [m]
                                  #     (of the bell function!)

db_NBUC   = 360.                  # NBUC thickness from core to bottom [m]
                                  #     (of the bell function!)

# BC-FM parameters =============================================================

BC_origin  = -14                  # BC Origin latitude

z0_BC      = 0.                   # BC core depth [m]

v0_BC_S    = -0.2                 # outgoing BC core velocity [m/s]

v0_BC_N    = 0                    # BC origin core velocity [m/s]

d_BC       = 150.                 # BC thickness [m] (of the bell function!)


# CLASSES ######################################################################


class ClassName(object):
    """docstring for ClassName"""
    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg
        

# FUNCTIONS ####################################################################

def transp(lon, z, v):
    dx = np.diff(lon, axis=1) * 111 * 1000  # valid only for low latitudes!!!
    aux  = dx[:,0]; aux.shape = (np.size(aux), 1)
    dx = np.concatenate( (dx, aux), axis=1) 

    dz = np.diff(z, axis=0)
    aux  = dz[0,:]; aux.shape = (1, np.size(aux))
    dz = np.concatenate( (dz, aux), axis=0)

    transp = np.abs(dx) * np.abs(dz) * v; transp = transp.sum()
    transp = transp / 1e6

    return transp
	



# ----------- ==================================================== ----------- #
#                             SCRIPT STARTS HERE                               #  
# ----------- ==================================================== ----------- #

# loading roms grid to get limits and topography
print ' '
print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
print ' '
# I need a bigger grid to get the isobath
grdfile  = nc.Dataset(largegridfile)

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h      = grdfile.variables['h'][:]

# getting an isobath
plt.figure(); con = plt.contour(lonr, latr, h, levels=[100] )
col = con.collections[0]; paths = col.get_paths()
path0 = paths[0]; isob = path0.vertices; plt.close('all')
# limiting isobath within model domain
f = np.where( (isob[:,1] >= lims[2]) & (isob[:,1] <= lims[3]) )
isob = isob[f[0],:]

# smoothing isobath
isob[:,0] = smooth(isob[:,0], window_len=201, window='hanning')
isob[:,1] = smooth(isob[:,1], window_len=101, window='hanning')

# now I load original small grid
grdfile  = nc.Dataset(smallgridfile)

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h      = grdfile.variables['h'][:]

# creating adimensional pairs of the parameters
dr    = dr / D                       # adimensional horizontal resolution  
r  = D * np.arange(0, 1, dr)         # adimensional horizontal transects
r0 = (1.0/D) * D                     # defining transect center
d  = ( delta / 111.0 ) / D           # normalized jet width


# ----------- ==================================================== ----------- #
#                             NBUC FEATURE MODEL                               #  
# ----------- ==================================================== ----------- #

# defining domain, buffering variables
looprange = range(0, len(isob), int(ln)) # scanning the isobath
li = np.size(looprange)
lj = np.size(r)
X  = np.zeros( [li , lj])
Y = X.copy()
Z  = np.arange(-zmax, 0.0+dz, dz)
lk = np.size(Z)

U = np.zeros( [ lk , li , lj ] )
V = U.copy(); VS = U.copy()
v = np.zeros( [ lk , lj ] )

# defining velocity-axis vertical structure v0 = v0(y,z)

# Y-dependency:  

# Jet core depth will increase as a linear function from south to north
z0max = np.linspace(-z0_NBUC_S, -z0_NBUC_N, li)

# Core velocity will also increase as a linear function from south to north
v0max = np.linspace(v0_NBUC_S, v0_NBUC_N, li) 

v0 = np.zeros( [lk, li] )

# Z-dependency:

# NBUC core in Z0 m decaying in a gauss curve until 0 m/s at bottom
#    this will also by Y-dependent, to allow increasing of the jet thickness
d1      = ds_NBUC / dz # gaussian vertical upper width, normalized by dz

# another gaussian will be adopted to decay velocities up to surface 
d2      = np.linspace(db_NBUC/dz, db_NBUC/dz, li) # gaussian lower width, 
                                                  # normalized by dz

# starting the loop to create the slope-following NBUC-FM
print ' '
print '======== CREATING SLOPE-FOLLOWING NBUC-FM ========'
print ' '
i = -1 # initializing index counter
for c in looprange:
    print '    Transect ' + str(i+1) + ' / ' + str(np.size(looprange))
    i = i + 1
    x0 = isob[c:c+6, 0] 
    y0 = isob[c:c+6, 1]
    tgr, b = np.polyfit(x0, y0, 1) # finding isobath-tangent straight line
    x0 = x0[0]
    y0 = y0[0]
    tgr = -1.0 / tgr
    b = y0 - tgr * x0 # finding normal straight line
    if tgr >= tgmax: # enforcing maximum angle
       tgr = tgmax
    elif tgr <= -tgmax:
       tgr = -tgmax
    
    # assembling jet-core vertical structure
    
    # upper part
    z1 = Z[z0max[i]/dz:]
    v01 = v0max[i] * np.exp( (-1)* (z1 - (z0max[i]/dz))**2 / (2*d1**2) )
    v0[z0max[i]/dz:, i] = v01

    # lower part
    z2  = Z[:z0max[i]/dz]
    v02 = v0max[i] * np.exp( (-1)* (z2 - (z0max[i]/dz))**2 / (2*d2[i]**2) )
    v0[:z0max[i]/dz, i] = v02
    
    # writing NBUC-FM to transect based on cross-slope structure
    
    for k in range( 0, lk):
        v[k, :] = v0[k, i] * np.exp( (-1)* ( ( r-r0 )**2 / ( 2*d**2 ) ) )

    # georreferencing velocities and coordinates

    angr = np.arctan(tgr)   
    cosr, sinr  = np.cos(angr), np.sin(angr)
    X[i, :]     = r * cosr + x0
    Y[i, :]     = r * sinr + y0
    U[:, i, :]  = v * sinr * (-1)
    V[:, i, :]  = v * cosr
    VS[:, i, :] = v
    
    
# ----------- ==================================================== ----------- #
#                             BC FEATURE MODEL                                 #  
# ----------- ==================================================== ----------- #

# defining domain, buffering variables
U2 = np.zeros( [ lk , li , lj ] )
V2 = U2.copy()
VS2 = U2.copy()
v2 = np.zeros( [ lk , lj ] )

# defining velocity-axis vertical structure v0 = v0(y,z)

# Y-dependency:  

# Core velocity will also increase as a linear function from south to north
v0max = np.zeros(li)
lataux = Y[:,0]
fcb = np.where(lataux <= BC_origin) 
fcb = fcb[0]
v0max[fcb] = np.linspace(v0_BC_S, v0_BC_N, fcb.size) 

v0 = np.zeros( [lk, li] )

# Z-dependency:

# NBUC core in Z0 m decaying in a gauss curve until 0 m/s at bottom
#    this will also by Y-dependent, to allow increasing of the jet thickness
d1      = d_BC / dz # gaussian vertical upper width, normalized by dz

# starting the looping to create the slope-following NBUC-FM
print ' '
print '======== CREATING SLOPE-FOLLOWING BC-FM ========'
print ' '
i = -1 # initializing index counter
for c in looprange:
    print '    Transect ' + str(i+1) + ' / ' + str(np.size(looprange))
    i = i + 1
    x0 = isob[c:c+6, 0]
    y0 = isob[c:c+6, 1]
    tgr, b = np.polyfit(x0, y0, 1) # finding isobath-tangent straight line
    x0 = x0[0]; y0 = y0[0]
    tgr = -1.0 / tgr
    b = y0 - tgr * x0 # finding normal straight line
    if tgr >= tgmax: # enforcing maximum angle
       tgr = tgmax
    elif tgr <= -tgmax:
       tgr = -tgmax
    
    # assembling vertical jet-core structure
    v0[:,i] = v0max[i] * np.exp( (-1)* (Z - (z0_BC/dz))**2 / (2*d1**2) )
    
    # writing NBUC-FM to transect based on cross-slope structure
    
    for k in range( 0, lk):
        v2[k, :] = v0[k, i] * np.exp( (-1)* ( ( r-r0 )**2 / ( 2*d**2 ) ) )

    # georreferencing velocities and coordinates

    angr = np.arctan(tgr)   
    cosr, sinr   = np.cos(angr), np.sin(angr)
    X[i, :]      = r * cosr + x0
    Y[i, :]      = r * sinr + y0
    U2[:, i, :]  = v2 * sinr * (-1)
    V2[:, i, :]  = v2 * cosr
    VS2[:, i, :] = v2
    
    
    
# Gathering both Feature Models
U = U + U2; V = V + V2; VS = VS + VS2
 

# ----------- ==================================================== ----------- #
#                             SOME PLOTTING                                    #  
# ----------- ==================================================== ----------- #

# loading map projections to make the plots
m = Basemap(projection='merc',
     llcrnrlat = latr.min(), urcrnrlat = latr.max(),
     llcrnrlon = lonr.min(), urcrnrlon = lonr.max(),
     lat_ts=0, resolution='l')
     
zb = np.ma.masked_where(h >= 100, h)
xb,yb = m(isob[:,0], isob[:,1])
lonm, latm = m(lonr, latr)   
Xm, Ym = m(X, Y)  

plt.figure(facecolor='w', figsize=(14,8))
# horizontal slices
p1=plt.subplot(1,2,1)
m.contourf(Xm, Ym, V[-1,...]*100, np.arange(-10, 10+0.3,0.3), 
           cmap=plt.cm.RdBu_r, extend='both')
plt.colorbar(shrink=0.9, extend='both')
m.contourf(lonm, latm, zb, colors=('0.7'), alpha=1)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,5], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,5], zorder=7)
tx, ty = m(lonr.min()+1, latr.max()-2)
txt = 'Velocities\n @ Surface\n [cm s$^{-1}$]'
p1.text(tx,ty,txt,color='k',fontsize=11,fontweight='bold')
bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = m(-33,-22)
plt.text(ax, ay, "a", ha="center", va="center", size=12, bbox=bbox_props)

p2=plt.subplot(1,2,2)
m.contourf(Xm, Ym, V[-400/dz,...]*100, np.arange(-50, 50+2,2), 
           cmap=plt.cm.RdBu_r, extend='both')
plt.colorbar(shrink=0.9, extend='both')
m.contourf(lonm, latm, zb, colors=('0.7'), alpha=1)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,5], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,5], zorder=7)
tx, ty = m(lonr.min()+1, latr.max()-2)
txt = 'Velocities\n @ 400m\n [cm s$^{-1}$]'
p2.text(tx,ty,txt,color='k',fontsize=11,fontweight='bold')
bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = m(-33,-22)
plt.text(ax, ay, "b", ha="center", va="center", size=12, bbox=bbox_props)

# vertical slice
pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(V[:,0,:])*100,np.arange(-20,20+3,3), 
                  cmap=plt.cm.RdBu, alpha=1)
plt.colorbar(shrink=0.9, extend='both')


# ----------- ==================================================== ----------- #
#                COMPUTING GEOSTROPHICALLY BALANCED MASS FIELD                 #  
# ----------- ==================================================== ----------- #

print ' '
print '======== COMPUTING GEOSTROPHICALLY BALANCED MASS FIELD ========'
print ' '
# integration the thermal wind equation:
# rho(x,z) = rho0(z) - rho_bar.f/g * int_0^L{ dv/dz dx}

# obtaining rho0 and rho_bar from WOA2009:
MeanDens = sp.loadmat(climatology_profile)
rho0  = MeanDens['potdens'][:].ravel()
rho0 = rho0[::-1]
zrho  = MeanDens['z'][:].ravel()
zrho = zrho[::-1]

salt0 = MeanDens['salt'][:].ravel()
salt0 = salt0[::-1] 

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



# ----------- ==================================================== ----------- #
#                                 SOME PLOTTING                                #  
# ----------- ==================================================== ----------- # 

# zonal slice
pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(RHO[:,0,:])-1000,
                  np.arange(23.2,27.8,0.2), cmap=plt.cm.Spectral_r, alpha=1)
plt.colorbar(shrink=0.9, extend='both')
plt.contour(X[0,:],Z, np.squeeze(RHO[:,0,:])-1000,
            np.arange(23.2,27.8,0.2), colors='0.4') 

# meridional slice 
plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
plt.contourf(Y[:,-1],Z, np.squeeze(RHO[:,:,-1])-1000,np.arange(23.2,27.8,0.2), 
             cmap=plt.cm.Spectral_r, alpha=1)
plt.colorbar(shrink=0.9, extend='both')
plt.contour(Y[:,-1],Z, np.squeeze(RHO[:,:,-1])-1000,
            np.arange(23.2,27.8,0.2), colors='0.4')


# ----------- ==================================================== ----------- #
#                   BUILDING THE REST OF THE VELOCITY FIELD                    #  
# ----------- ==================================================== ----------- #

print ' '
print '======== BUILDING THE REST OF THE VELOCITY FIELD ========'
print ' '

# clearing memory
del h, VS, MeanDens, U2, V2, VS2
del c, col, con, cosr, d, d1, d2, dvdz, fcb, grdfile

# extrapolate to the rest of the domain 
rpt = 20
XAUX = np.zeros([li, rpt]); YAUX = XAUX * 0

print ' '
print '    Extrapolating values to the east'
print ' '
for i in range(0,li):
    lastx = X[i,-1]; xaux = np.linspace(lastx+0.25, lonr.max()+2, rpt)
    XAUX[i,:] = xaux
    lasty = Y[i,-1]; yaux = np.linspace(lasty, lasty, rpt)
    YAUX[i,:] = yaux

# coordinates:
X  = np.hstack((X, XAUX))
Y  = np.hstack((Y, YAUX))

# velocity:

# computing geostrophic velocity at the EASTERN boundary, to get SEC structure

temp = np.squeeze(TEMP[...,-1])
salt = np.squeeze(SALT[...,-1])
y    = Y[:,-1] 
y, z = np.meshgrid(y,Z); z = -z

gp   = sw.gpan(salt, temp, z)
gp   = (gp - gp[-1,:]) * -1 # relative velocities in respect to the bottom
 
dgp  = np.array(np.gradient(gp))
dgp  = np.squeeze(dgp[1,:,:])

dy   = np.array(np.gradient(y))
dy   = np.squeeze(dy[1,:,:]) * 111000

dgpdy = dgp / dy 
usec  = -dgpdy / f0

usec = usec*0.95 # getting the right transport
f = np.where(usec > 0); usec[f] = 0

Tsec = transp(y, z, usec)

# checking SEC's vertical structure and transport
plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
plt.contourf(y, -z, usec*100, np.arange(-10, 10+1, 1), cmap=plt.cm.RdBu);
plt.colorbar(shrink=0.9, extend='both')
plt.axis([latr.min(), latr.max(), -1250, 50])
plt.title('East')
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
plt.text(latr.mean(), -400, r'SEC Transport = 17 $Sv$',
         ha="center", va="center", size=16, bbox=bbox_props)
         

# keep going with final fields construction - extrapolations
uaux = usec.repeat(rpt, axis=1)
uaux.shape = (lk, li, rpt)

# decaying uaux to zero close to the western boundary
dcj = 6 # number of grid points to do the linear decay

for k in range(0, lk-1):
    for i in range(0, li):
        utmp = np.linspace(0, usec[k, i], dcj)
        uaux[k, i, :dcj] = utmp
          

U  = np.concatenate((U, uaux), axis=2)
V  = np.concatenate((V, uaux*0), axis=2)

del uaux

# salt and temp:
lastS, lastT, lastR = SALT[:,:,-1], TEMP[:,:,-1], RHO[:,:,-1]

saux = lastS.repeat(rpt, axis=1)
taux = lastT.repeat(rpt, axis=1)
raux = lastR.repeat(rpt, axis=1)

saux.shape = (lk, li, rpt) 
taux.shape = (lk, li, rpt)
raux.shape = (lk, li, rpt)

SALT  = np.concatenate((SALT, saux), axis=2)
TEMP  = np.concatenate((TEMP, taux), axis=2)
RHO   = np.concatenate((RHO,  raux), axis=2)

lk, li, lj = U.shape

rpt = 10
XAUX = np.zeros([li, rpt]); YAUX = XAUX * 0
print ' '
print '    Extrapolating values to the west'
print ' '
for i in range(0,li):
    firstx = X[i,0]
    xaux = np.linspace(lonr.min() - 2, firstx-0.25, rpt)
    XAUX[i,:] = xaux
    
    firsty = Y[i,0] 
    yaux = np.linspace(firsty, firsty, rpt)
    YAUX[i,:] = yaux

X  = np.hstack( (XAUX, X) )
Y  = np.hstack( (YAUX, Y) )

firstu = U[:,:,0]
firstu = firstu*0 
uaux = firstu.repeat(rpt, axis=1)
uaux.shape = (lk, li, rpt)
U  = np.concatenate((uaux, U), axis=2)
V  = np.concatenate((uaux, V), axis=2)

firstS, firstT, firstR = SALT[:,:,0], TEMP[:,:,0], RHO[:,:,0]  

saux = firstS.repeat(rpt, axis=1)
taux = firstT.repeat(rpt, axis=1)
raux = firstR.repeat(rpt, axis=1)

saux.shape, taux.shape, raux.shape = (lk, li, rpt), (lk, li, rpt), (lk, li, rpt)

SALT  = np.concatenate((saux, SALT), axis=2)
TEMP  = np.concatenate((taux, TEMP), axis=2)
RHO   = np.concatenate((raux,  RHO), axis=2)


# ----------- ==================================================== ----------- #
#                                 SOME PLOTTING                                #  
# ----------- ==================================================== ----------- #

Xm, Ym = m(X, Y)  
plt.figure(facecolor='w', figsize=(14,8))

p1 = plt.subplot(1,2,1)
m.contourf(Xm, Ym, RHO[-1,...]-1000, np.arange(23.5,24,0.03), 
           cmap=plt.cm.Spectral_r, alpha=1)
plt.colorbar(shrink=0.8, extend='both')
m.contour(Xm, Ym, RHO[-1,...]-1000, np.arange(23.5,24,0.03), colors='0.4')
m.contourf(lonm, latm, zb, colors=('0.9'), alpha=0.5)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
tx, ty = m(lonr.min()+1, latr.max()-2)
text = '$\sigma_\Theta$\n @ Surface\n [$kg.m^{-3}$]'
p1.text(tx,ty,text,color='k',fontsize=11,fontweight='bold')

p2 = plt.subplot(1,2,2)
m.contourf(Xm, Ym, RHO[-400/dz,...]-1000, np.arange(26.28,27.06,0.03), 
           cmap=plt.cm.Spectral_r, alpha=1)
plt.colorbar(shrink=0.8, extend='both')
m.contour(Xm, Ym, RHO[-400/dz,...]-1000, 
          np.arange(26.28,27.02,0.03), colors='0.4')
m.contourf(lonm, latm, zb, colors=('0.9'), alpha=0.5)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
tx, ty = m(lonr.min()+1, latr.max()-2)
text = '$\sigma_\Theta$\n @ 400 m\n [$kg.m^{-3}$]'
p2.text(tx,ty,text,color='k',fontsize=11,fontweight='bold')



# ----------- ==================================================== ----------- #
#                          PREPARING FIELDS FOR ROMS                           #  
# ----------- ==================================================== ----------- #

print ' '
print '======== PREPARING FIELDS FOR ROMS ========'
print ' '

# clearing memory
del lastS, lastT, lastx, lasty, firstS, firstT, firstu, firstx, firsty
del lataux, looprange, r, rho0, rhoaux, taux, saux, tmp, v, v0, v01, v02, v0max
del XAUX, YAUX, aux, uaux

d = 20
TEMP, SALT, Z, U, V = TEMP[::d,...], SALT[::d,...], Z[::d], U[::d,...], V[::d,...]
lk, li, lj = U.shape

# flipping arrays to keep depths positive
TEMP, SALT = TEMP.reshape(lk, li*lj), SALT.reshape(lk, li*lj)
U, V = U.reshape(lk, li*lj), V.reshape(lk, li*lj)

TEMP, SALT = np.flipud(TEMP), np.flipud(SALT)
U, V = np.flipud(U), np.flipud(V)

Z = np.flipud(Z); Z = -Z; Z = np.squeeze(Z)

TEMP.shape = (lk,li,lj)
SALT.shape = (lk,li,lj)
U.shape = (lk,li,lj)
V.shape = (lk,li,lj)

# creating depth-averaged velocities
UBAR = U.mean(axis=0)
VBAR = V.mean(axis=0)

# creating SSH
TEMP2 = TEMP.reshape(lk, li*lj)
SALT2 = SALT.reshape(lk, li*lj)
Z2    = Z.reshape(lk,1)
gpan  = sw.gpan(SALT2, TEMP2, Z2)
gpan  = (gpan - gpan[-1,:]) # to reference in the bottom
SSH   = gpan / g 
SSH = SSH[0,:]
SSH = SSH - SSH.mean()
SSH.shape = (li,lj)


# ----------- ==================================================== ----------- #
#                                 SOME PLOTTING                                #  
# ----------- ==================================================== ----------- #


plt.figure(facecolor='w', figsize=(14,8))

p1 = plt.subplot(1,2,1)
m.contourf(Xm, Ym, SSH*100, 30, cmap=plt.cm.RdBu_r, extend='both')
plt.colorbar(shrink=0.8, extend='both')
#m.contour(Xm, Ym, SSH*-100,30, colors='0.4')
m.contourf(lonm, latm, zb, colors=('0.9'), alpha=0.5)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
tx, ty = m(lonr.min()+1, latr.max()-1)
text = 'SSH [$cm$]'
p1.text(tx,ty,text,color='k',fontsize=11,fontweight='bold')
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)

p2 = plt.subplot(1,2,2)
m.contourf(Xm, Ym, UBAR*100, np.arange(-2, 2+0.1,0.1), cmap=plt.cm.RdBu_r, 
           extend='both')
plt.colorbar(shrink=0.8, extend='both')
m.contourf(lonm, latm, zb, colors=('0.9'), alpha=0.5)
m.plot(xb, yb,'k--');
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latr.min(), latr.max(), 2),
   labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
   labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
tx, ty = m(lonr.min()+0.6, latr.max()-2)
text = 'Depth-averaged\n zonal velocities\n @ 400 m\n [$cm.s^{-1}$]'
p2.text(tx,ty,text,color='k',fontsize=11,fontweight='bold')



print ' '
print '======== PLOTTING ========'
print ' '
plt.show()



# ----------- ==================================================== ----------- #
#                             SCRIPT ENDS HERE                                 #  
# ----------- ==================================================== ----------- #








