#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Building coupled BC-NBUC Feature Model
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Apr, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from roms_setup import near 
import scipy.io as sp
from mpl_toolkits.basemap import Basemap 
import seawater.csiro as sw
import netCDF4 as nc
from cookb_signalsmooth import smooth

# FUNCTIONS #################################################################

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
	


# Input Parameters ##################################################################

# Common System Parameters =========================================================

D         = 3.                    # length of the transect or total length of the 
                                  # modeled jet [degrees]
                                     
tgmax     = 0.15                   # tan of the maximum angle that the transect is
                                  # allowed to have to a parallel
                                     
dr        = 0.1                  # horizontal spacing for the feature model [degree]
                                     
ln        = 30.                   # jump on the isobath to define the # of transects

zmax      = 5000.                 # maximum depth for the system

dz        = 1.                    # vertical resolution of the transects [m]
                                  # DO NOT CHANGE THAT!!!! If changed, 
                                  # require changes in the code

delta     = 100.                  # jets width [km]

# NBUC parameters ==================================================================

z0_NBUC_S = 500.                  # incoming NBUC core depth [m]

z0_NBUC_N = 200.                  # outgoing NBUC core depth [m]

v0_NBUC_S = 0.0002                   # incoming NBUC core velocity [m/s]

v0_NBUC_N = 0.0005                   # outgoing NBUC core velocity [m/s]

ds_NBUC   = 100.                  # NBUC thickness from core to surface [m]

db_NBUC   = 360.                  # NBUC thickness from core to bottom [m]

# BC parameters ===================================================================

BC_origin  = -14                  # BC Origin latitude

z0_BC      = 0.                   # BC core depth [m]

v0_BC_S    = -0.2              # outgoing BC core velocity [m/s]

v0_BC_N    = 0                    # BC origin core velocity [m/s]

d_BC       = 150.                 # BC thickness [m] 

#####################################################################################

# analytical integration of the functions to check on the transport

# NBUC at North
## Analytical Expression of the Transport: T = v0 . d ( ds + db )

#plt.figure(10, figsize=(14,8), facecolor='w')
## v0, d, ds fixed:
#db = np.arange(150, 500, 10)
#T  =  ( v0_NBUC_N * d*1e3 * ( ds_NBUC + db ) ) * 1e-6
#p1 = plt.subplot(221); plt.plot(db, T, 'g', linewidth=2)
#plt.plot( (150,500) , (23,23) , 'k', alpha='0.2', linewidth=10); grid()
#plt.title('$v_0$ = '+ str(v0_NBUC_N) +'m/s, $\delta$ = '+str(int(d)) +' km, $\delta_s$ = '+str(int(ds_NBUC)) +' m')
#plt.xlabel('$\delta_b$ [m]'); plt.ylabel('NBUC Transport [Sv]'); p1.set_ylim(15,30)

#ds = np.arange(20, 300, 10)
#T  =  ( v0_NBUC_N * d*1e3 * ( ds + db_NBUC ) ) * 1e-6
#p2 = plt.subplot(222); plt.plot(ds, T, 'g', linewidth=2)
#plt.plot( (20,300) , (23,23) , 'k', alpha='0.2', linewidth=10); grid()
#plt.title('$v_0$ = '+ str(v0_NBUC_N) +'m/s, $\delta$ = '+str(int(d)) +' km, $\delta_b$ = '+str(int(db_NBUC)) +' m')
#plt.xlabel('$\delta_s$ [m]'); p2.set_ylim(15,30)

#dd = np.arange(50, 150, 10)
#T  =  ( v0_NBUC_N * dd*1e3 * ( ds_NBUC + db_NBUC ) ) * 1e-6
#p3 = plt.subplot(223, position=(0.125, 0.1, 0.35, 0.3)); plt.plot(dd, T, 'g', linewidth=2)
#plt.plot( (50,150) , (23,23) , 'k', alpha='0.2', linewidth=10); grid()
#plt.title('$v_0$ = '+ str(v0_NBUC_N) +'m/s, $\delta_s$ = '+str(int(ds_NBUC)) +' m, $\delta_b$ = '+str(int(db_NBUC)) +' m')
#plt.xlabel('$\delta$ [km]');  plt.ylabel('NBUC Transport [Sv]'); p3.set_ylim(15,30); 

#v0 = np.arange(0.3, 0.8, 0.03)
#T  =  ( v0 * d*1e3 * ( ds_NBUC + db_NBUC ) ) * 1e-6
#p4 = plt.subplot(224, position=(0.55, 0.1, 0.35, 0.3)); plt.plot(v0, T, 'g', linewidth=2)
#plt.plot( (0.3,0.8) , (23,23) , 'k', alpha='0.2', linewidth=10); grid()
#plt.title('$\delta$ = '+ str(int(d)) +'km, $\delta_s$ = '+str(int(ds_NBUC)) +' m, $\delta_b$ = '+str(int(db_NBUC)) +' m')
#plt.xlabel('$v_0$ [m/s]'); p4.set_ylim(15,30)



# ======================================================
# CREATING ISOBATH-FOLLOWING NBUC FEATURE MODEL:
# ======================================================

# ======================================================
# loading roms grid to get limits and topography
print ' '
print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
print ' '
# I need a bigger grid to get the isobath
grdfile  = nc.Dataset('/home/rsoutelino/myroms/phd_run/phd1_grd.nc')

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h      = grdfile.variables['h'][:]

# getting an isobath
plt.figure(); con = plt.contour(lonr, latr, h, levels=[100] )
col = con.collections[0]; paths = col.get_paths()
path0 = paths[0]; isob = path0.vertices; plt.close('all')
# limiting isobath within model domain
f = np.where( (isob[:,1] >= -24) & (isob[:,1] <= -8) )
isob = isob[f[0],:]

# smoothing isobath
isob[:,0] = smooth(isob[:,0],window_len=201,window='hanning')
isob[:,1] = smooth(isob[:,1],window_len=101,window='hanning')

# now I load original small grid
grdfile  = nc.Dataset('/home/rsoutelino/myroms/phd_run/phd8_grd.nc')

# assigning some variables from grid file
lonr   = grdfile.variables['lon_rho'][:]
latr   = grdfile.variables['lat_rho'][:]
h      = grdfile.variables['h'][:]

# creating adimensional pairs of the parameters
dr    = dr / D                       # adimensional horizontal resolution  
r  = D * np.arange(0, 1, dr)         # adimensional horizontal transects
r0 = (1.0/D) * D                     # defining transect center
d  = ( delta / 111.0 ) / D           # normalized jet width

### NBUC FEATURE MODEL #####################################################
# ======================================================
# defining domain, buffering variables
looprange = range(0, len(isob), int(ln)) # scanning the isobath
li = np.size(looprange); lj = np.size(r)
X  = np.zeros( [li , lj]); Y = X.copy()
Z  = np.arange(-zmax, 0.0+dz, dz); lk = np.size(Z)

U = np.zeros( [ lk , li , lj ] )
V = U.copy(); VS = U.copy()
v = np.zeros( [ lk , lj ] )

# ======================================================
# defining velocity-axis vertical structure v0 = v0(y,z)

# Y-dependance:  

# Jet core depth will increase as a linear function from south to north
z0max = np.linspace(-z0_NBUC_S, -z0_NBUC_N, li)

# Core velocity will also increase as a linear function from south to north
v0max = np.linspace(v0_NBUC_S, v0_NBUC_N, li) 

v0 = np.zeros( [lk, li] )

# Z-dependance:
# NBUC core in Z0 m decaying in a gauss curve until 0 m/s at bottom
#    this will also by Y-dependant, to allow increasing of the jet thickness
d1      = ds_NBUC / dz # gaussian vertical upper width, normalized by dz

# another gaussian will be adopted to decay velocities to surface 
d2      = np.linspace(db_NBUC/dz, db_NBUC/dz, li) # gaussian lower width, normalized by dz

# starting the looping to create the slope-following NBUC-FM
print ' '
print '======== CREATING SLOPE-FOLLOWING NBUC-FM ========'
print ' '
i = -1 # initializing index counter
for c in looprange:
    print '    Transect ' + str(i+1) + ' / ' + str(np.size(looprange))
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
    
    
### BC FEATURE MODEL #####################################################
# ======================================================
# defining domain, buffering variables
U2 = np.zeros( [ lk , li , lj ] )
V2 = U2.copy(); VS2 = U2.copy()
v2 = np.zeros( [ lk , lj ] )

# ======================================================
# defining velocity-axis vertical structure v0 = v0(y,z)

# Y-dependance:  

# Core velocity will also increase as a linear function from south to north
v0max = np.zeros(li)
lataux = Y[:,0]
fcb = np.where(lataux <= BC_origin); fcb = fcb[0]
v0max[fcb] = np.linspace(v0_BC_S, v0_BC_N, fcb.size) 

v0 = np.zeros( [lk, li] )

# Z-dependance:
# NBUC core in Z0 m decaying in a gauss curve until 0 m/s at bottom
#    this will also by Y-dependant, to allow increasing of the jet thickness
d1      = d_BC / dz # gaussian vertical upper width, normalized by dz

# starting the looping to create the slope-following NBUC-FM
print ' '
print '======== CREATING SLOPE-FOLLOWING BC-FM ========'
print ' '
i = -1 # initializing index counter
for c in looprange:
    print '    Transect ' + str(i+1) + ' / ' + str(np.size(looprange))
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
    v0[:,i] = v0max[i] * np.exp( (-1)* (Z - (z0_BC/dz))**2 / (2*d1**2) )
    
    # ==========================================================
    # writing NBUC-FM to transect based on cross-slope structure
    
    for k in range( 0, lk):
        v2[k, :] = v0[k, i] * np.exp( (-1)* ( ( r-r0 )**2 / ( 2*d**2 ) ) )

    # georeferencing velocities and coordinates
    angr = np.arctan(tgr)   
    cosr, sinr   = np.cos(angr), np.sin(angr)
    X[i, :]      = r * cosr + x0
    Y[i, :]      = r * sinr + y0
    U2[:, i, :]  = v2 * sinr * (-1)
    V2[:, i, :]  = v2 * cosr
    VS2[:, i, :] = v2
    
    
    
# Gathering both Feature Models
U = U + U2; V = V + V2; VS = VS + VS2
 

# some plotting and transport computation

plt.figure() 
plt.plot(isob[:,0], isob[:,1]);  plt.axis('equal')
plt.pcolormesh(X, Y, V[-1,...], vmin=-0.20, vmax=0.20, cmap=plt.cm.RdBu)
plt.grid(); plt.colorbar(); plt.title('V-Vel @ Surface')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure() 
plt.plot(isob[:,0], isob[:,1]);  plt.axis('equal')
plt.pcolormesh(X, Y, V[-400/dz,...], vmin=-0.50, vmax=0.50, cmap=plt.cm.RdBu)
plt.grid(); plt.colorbar(); plt.title('V-Vel @ 400 m')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure()
plt.contourf(X[0,:],Z, np.squeeze(V[:,0,:]),np.arange(-0.2,0.2+0.03,0.03), cmap=plt.cm.RdBu, alpha=0.5)
plt.colorbar(); plt.title('V-Vel @ South Boundary')
plt.xlabel('Longitude'); plt.ylabel('Depth')
xx,zz = np.meshgrid(X[0,:],Z)
Tcb   = transp(xx[-200:,:], zz[-200:,:], np.squeeze(V[-200:,0,:]))
Tnbuc = transp(xx[:-200,:], zz[:-200,:], np.squeeze(V[:-200,0,:]))
plt.text(-42.5, -100,'BC: '+str(np.round(Tcb))+' Sv',color='k',fontsize=12,fontweight='bold')
plt.text(-42.5, -700,'NBUC: '+str(np.round(Tnbuc))+' Sv',color='k',fontsize=12,fontweight='bold')

###############################################

plt.figure()
plt.contourf(X[-1,:],Z, np.squeeze(V[:,-1,:]),np.arange(-0.5,0.5+0.05,0.05), cmap=plt.cm.RdBu, alpha=0.5)
fwhm_s = 2*np.sqrt(2*np.log1p(2)) * ds_NBUC; zd = (-z0_NBUC_N, -z0_NBUC_N + fwhm_s/2 )
plt.plot((-34.13,-34.13),zd,'k', linewidth=5, alpha=0.4)
plt.text(-34, -90, '$\delta_s$', fontsize=16, fontweight='bold')
plt.text(-33.92, -90, ' = '+ str(int(ds_NBUC)) +'m')
fwhm_b = 2*np.sqrt(2*np.log1p(2)) * db_NBUC; zd = (-z0_NBUC_N, -z0_NBUC_N - fwhm_b/2 )
plt.plot((-34.13,-34.13),zd,'b', linewidth=5, alpha=0.4)
plt.text(-34, -500, '$\delta_b$', fontsize=16, fontweight='bold', color='b')
plt.text(-33.92, -500, ' = '+ str(int(db_NBUC)) +'m', fontsize=12, color='b')
fwhm = delta/111; xd = (-34.13 + fwhm/2, -34.13 - fwhm/2)
plt.plot(xd,(-200,-200),'g', linewidth=5, alpha=0.4)
plt.text(-33.8, -250, '$\delta$', fontsize=16, fontweight='bold', color='g')
plt.text(-33.72, -250, ' = '+ str(int(delta)) +'km', fontsize=12, color='g')
plt.text(-34.2, -220, '$v_0$', fontsize=20, fontweight='bold', color='k')
plt.text(-32.9, -250, '$v_0$', fontsize=20, fontweight='bold', color='k')
plt.text(-32.78, -250, ' = '+ str(v0_NBUC_N) +' m/s', fontsize=12, color='k')
plt.colorbar(); plt.title('V-Vel @ North Boundary')
plt.xlabel('Longitude'); plt.ylabel('Depth')
xx,zz = np.meshgrid(X[0,:],Z)
Tnbuc = transp(xx, zz, np.squeeze(V[:,-1,:]))
plt.text(-33, -200,'NBUC: '+str(np.round(Tnbuc))+' Sv',color='k',fontsize=12,fontweight='bold')
plt.axis([-35, -32, -5000, 10])
plt.show() 
    
   





    
    
   
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


plt.figure() 
plt.plot(isob[:,0], isob[:,1]);  plt.axis('equal')
plt.pcolormesh(X, Y, TEMP[-1,...],cmap=plt.cm.Spectral_r); 
plt.grid(); plt.colorbar(); plt.title('Temperature @ Surface')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure() 
plt.plot(isob[:,0], isob[:,1]);  plt.axis('equal')
plt.pcolormesh(X, Y, TEMP[-400,...], cmap=plt.cm.Spectral_r)
plt.grid(); plt.colorbar(); plt.title('Temperature @ 400 m')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure()
plt.contourf(X[0,:],Z, np.squeeze(TEMP[:,0,:]),30, cmap=plt.cm.Spectral_r); plt.colorbar(); 
plt.contour(X[0,:],Z, np.squeeze(V[:,0,:]),10, colors='0.5')
plt.title('Temperature @ South Boundary')
plt.xlabel('Longitude'); plt.ylabel('Depth')

plt.figure()
plt.contourf(X[-1,:],Z, np.squeeze(TEMP[:,-1,:]),30, cmap=plt.cm.Spectral_r); plt.colorbar()
plt.contour(X[-1,:],Z, np.squeeze(V[:,-1,:]),10, colors='0.5')
plt.title('Temperature @ North Boundary')
plt.xlabel('Longitude'); plt.ylabel('Depth')
plt.show()

plt.figure()
plt.contourf(Y[:,-1],Z, np.squeeze(TEMP[:,:,-1]),30, cmap=plt.cm.Spectral_r); plt.colorbar()
plt.title('Temperature @ East Boundary')
plt.xlabel('Latitude'); plt.ylabel('Depth')
plt.show()
 
# clearing memory
del h, VS, MeanDens, RHO, U2, V2, VS2
del c, col, con, cosr, d, d1, d2, dvdz, fcb, grdfile


# extrapolate to the rest of the domain 
rpt = 20
XAUX = np.zeros([li, rpt]); YAUX = XAUX * 0

print ' '
print ' Extrapolating values to the east'
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

# computing geostrophic velocity to EASTERN boundary, to check SEC structure
temp = np.squeeze(TEMP[...,-1])
salt = np.squeeze(SALT[...,-1])
y    = Y[:,-1] 
y, z = np.meshgrid(y,Z); z = -z

gp   = sw.gpan(salt, temp, z)
gp   = (gp - gp[-1,:]) * -1 # to reference in the bottom
 
dgp  = np.array(np.gradient(gp))
dgp  = np.squeeze(dgp[1,:,:])

dy   = np.array(np.gradient(y))
dy   = np.squeeze(dy[1,:,:]) * 111000

dgpdy = dgp / dy 

usec  = -dgpdy / f0

# getting the right transport
usec = usec*150
f = np.where(usec > 0); usec[f] = 0

Tsec = transp(y, z, usec)

plt.figure(20, figsize=(10, 5), facecolor='w')
plt.contourf(y, -z, usec, np.arange(-0.1, 0.1+0.01, 0.01), cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.title('SEC velocities from BC-NBUC Feature Model')
plt.xlabel('Latitude'); plt.ylabel('Z[m]')
plt.axis([y.min(), y.max(), -5000, 0])
plt.text(-18, -450, str(np.round(Tsec))+' Sv')
plt.text(-23, -700, 'Should be -4 Sv to fulfil the BC-NBUC inbalance.')
plt.show()
 
uaux = usec.repeat(rpt, axis=1)
uaux.shape = (lk, li, rpt)

# decaying uaux to zero close to the western boundary
dcj = 3 # number of grid points to do the linear decay

for k in range(0, lk-1):
    for i in range(0, li):
        utmp = np.linspace(0, usec[k, i], dcj)
        uaux[k, i, :dcj] = utmp
          

U  = np.concatenate((U, uaux), axis=2)
V  = np.concatenate((V, uaux*0), axis=2)

del uaux

# salt and temp:
lastS = SALT[:,:,-1]; lastT = TEMP[:,:,-1]  
saux = lastS.repeat(rpt, axis=1); taux = lastT.repeat(rpt, axis=1)
saux.shape = (lk, li, rpt); taux.shape = (lk, li, rpt)
SALT  = np.concatenate((SALT, saux), axis=2)
TEMP  = np.concatenate((TEMP, taux), axis=2)


lk, li, lj = U.shape

rpt=10
XAUX = np.zeros([li, rpt]); YAUX = XAUX * 0
print ' '
print ' Extrapolating values to the west'
print ' '
for i in range(0,li):
    firstx = X[i,0]; xaux = np.linspace(lonr.min() - 2, firstx-0.25, rpt)
    XAUX[i,:] = xaux
    firsty = Y[i,0]; yaux = np.linspace(firsty, firsty, rpt)
    YAUX[i,:] = yaux

X  = np.hstack((XAUX, X))
Y  = np.hstack((YAUX, Y))

firstu = U[:,:,0]; firstu = firstu*0; 
uaux = firstu.repeat(rpt, axis=1)
uaux.shape = (lk, li, rpt)
U  = np.concatenate((uaux, U), axis=2)
V  = np.concatenate((uaux, V), axis=2)

firstS = SALT[:,:,0]; firstT = TEMP[:,:,0]  
saux = firstS.repeat(rpt, axis=1); taux = firstT.repeat(rpt, axis=1)
saux.shape = (lk, li, rpt); taux.shape = (lk, li, rpt)
SALT  = np.concatenate((saux, SALT), axis=2)
TEMP  = np.concatenate((taux, TEMP), axis=2)


del lastS, lastT, lastx, lasty, firstS, firstT, firstu, firstx, firsty
del lataux, looprange, r, rho0, rhoaux, taux, saux, tmp, v, v0, v01, v02, v0max
del v2, vaux, xx, yaux, z0max, z1, z2, zrho, zz, XAUX, YAUX, aux, lonr, latr, uaux

d = 20
TEMP = TEMP[::d,...]; SALT = SALT[::d,...]; Z = Z[::d]
U    =    U[::d,...];    V =    V[::d,...]
lk, li, lj = U.shape


# PREPARING FIELDS TO BUILD ROMS INITIAL FIELDS straight from FM:
# flipping arrays to keep depths positive
TEMP = TEMP.reshape(lk, li*lj); SALT = SALT.reshape(lk, li*lj)
U = U.reshape(lk, li*lj); V = V.reshape(lk, li*lj)
TEMP = np.flipud(TEMP); SALT = np.flipud(SALT); U = np.flipud(U); V = np.flipud(V);
Z = np.flipud(Z); Z = -Z; Z = np.squeeze(Z)
TEMP.shape = (lk,li,lj); SALT.shape = (lk,li,lj); U.shape = (lk,li,lj); V.shape = (lk,li,lj)

# creating depth-averaged velocities
UBAR = U.mean(axis=0)
VBAR = V.mean(axis=0)

# creating SSH
TEMP2 = TEMP.reshape(lk, li*lj)
SALT2 = SALT.reshape(lk, li*lj)
Z2    = Z.reshape(lk,1)
gpan  = sw.gpan(SALT2, TEMP2, Z2); del TEMP2, SALT2, Z2
gpan  = (gpan - gpan[-1,:]) # to reference in the bottom
SSH   = gpan / g; SSH = SSH[0,:]; SSH = SSH - SSH.mean()
SSH.shape = (li,lj)

matdict = {'lon':X, 'lat': Y, 'z': Z, 'temp':TEMP, 'salt':SALT, 'u':U, 'v':V, 'ubar':UBAR, 'vbar':VBAR, 'ssh':SSH}
sp.savemat('BC_FM.mat', matdict)

# CREATING A FLAT Tracers FM + M3 FM velocities
for k in range(0,lk-1):
    TEMP[k,...] = TEMP[k,...].mean()
    SALT[k,...] = SALT[k,...].mean()

matdict = {'lon':X, 'lat': Y, 'z': Z, 'temp':TEMP, 'salt':SALT, 'u':U, 'v':V, 'ubar':UBAR, 'vbar':VBAR, 'ssh':SSH}
sp.savemat('FLAT-BC-NBUC_FM.mat', matdict)

# CREATING A FLAT Tracers FM + FLAT M3 FM 
U = U*0; UBAR = UBAR*0
V = V*0; VBAR = VBAR*0; SSH = SSH*0

matdict = {'lon':X, 'lat': Y, 'z': Z, 'temp':TEMP, 'salt':SALT, 'u':U, 'v':V, 'ubar':UBAR, 'vbar':VBAR, 'ssh':SSH}
sp.savemat('FLAT_FM.mat', matdict)

# PREPARING FIELDS TO AOME
print ' '
print 'Please run fm2mod.py if you want to create a MODS-file' 
print ' '

















STOP
# comparing with OA fields

dataset = nc.Dataset('/home/rsoutelino/myroms/phd_run/init/hops_oa/work/bc-nbuc_fm.nc')
lon = dataset.variables['grid3'][:,:,0]
lat = dataset.variables['grid3'][:,:,1]
temp = dataset.variables['temp'][:]
salt = dataset.variables['salt'][:]
dynht = dataset.variables['dynht'][:]
z = dataset.variables['zout'][:,2]

dynht = dynht[0,...]
psi = (-1) * ( dynht / f0 )
gradpsi = np.gradient(psi)

u = (-1) * ( gradpsi[0] / ( np.diff(lon).mean() * 111000 ) )
v = gradpsi[1] / ( np.diff(lon).mean() * 111000 )

plt.figure() 
plt.contourf(lon, lat, dynht[...,0], 30, cmap=plt.cm.RdBu); colorbar()
plt.plot(isob[:,0], isob[:,1],'k', linewidth=2); axis('equal')
plt.axis([-40, -31, -23, -10])
plt.title('$\Delta \Phi$ @ Surface')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure() 
plt.contourf(lon, lat, dynht[...,20], 30, cmap=plt.cm.RdBu); plt.colorbar()
plt.plot(isob[:,0], isob[:,1],'k', linewidth=2); axis('equal')
plt.axis([-40, -31, -23, -10])
plt.title('$\Delta \Phi$ @ 400 m')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure() 
plt.quiver(lon, lat, u[...,0], v[...,0])
plt.plot(isob[:,0], isob[:,1],'k', linewidth=2);
plt.axis([-40, -31, -23, -10])
plt.title('Velocity @ Surface')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure() 
plt.quiver(lon, lat, u[...,20], v[...,20])
plt.plot(isob[:,0], isob[:,1],'k', linewidth=2);
plt.axis([-40, -31, -23, -10])
plt.title('Velocity @ 400 m')
plt.xlabel('Longitude'); plt.ylabel('Latitude')

plt.figure()
plt.contourf(X[0,:],Z, np.squeeze(V[:,0,:]),np.arange(-0.2,0.2+0.03,0.03), cmap=plt.cm.RdBu, alpha=0.5)
plt.contourf(lat[:,-1],z, np.squeeze(u[:,-60,:]).transpose(),30, cmap=plt.cm.RdBu); plt.colorbar()
plt.title('Temperature @ South Boundary')
plt.xlabel('Longitude'); plt.ylabel('Depth')
plt.show()















