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
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

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

zmax      = 1500.                 # maximum depth for the system

dz        = 1.                    # vertical resolution of the transects [m]
                                  # DO NOT CHANGE THAT!!!! If changed, 
                                  # require changes in the code

delta     = 100.                  # jets width [km]

# NBUC parameters ==================================================================

z0_NBUC_S = 500.                  # incoming NBUC core depth [m]

z0_NBUC_N = 200.                  # outgoing NBUC core depth [m]

v0_NBUC_S = 0.2                   # incoming NBUC core velocity [m/s]

v0_NBUC_N = 0.5                   # outgoing NBUC core velocity [m/s]

ds_NBUC   = 100.                  # NBUC thickness from core to surface [m]

db_NBUC   = 360.                  # NBUC thickness from core to bottom [m]

# BC parameters ===================================================================

BC_origin  = -14                  # BC Origin latitude

z0_BC      = 0.                   # BC core depth [m]

v0_BC_S    = -0.2                 # outgoing BC core velocity [m/s]

v0_BC_N    = 0                    # BC origin core velocity [m/s]

d_BC       = 150.                  # BC thickness [m] 

#####################################################################################

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
 
m = Basemap(projection='cyl',
     llcrnrlat = latr.min(), urcrnrlat = latr.max(),
     llcrnrlon = lonr.min(), urcrnrlon = lonr.max(),
     lat_ts=0, resolution='l')
     
zb = np.ma.masked_where(h >= 100, h)
xb,yb = m(isob[:,0], isob[:,1])
lonm, latm = m(lonr, latr)   
Xm, Ym = m(X, Y)  
     
#############################################################################################
plt.figure(facecolor='w', figsize=(14,8))

p1=plt.subplot(1,2,1)
m.contourf(Xm, Ym, V[-1,...]*100, np.arange(-10, 10+0.3,0.3), cmap=plt.cm.RdBu_r, extend='both')
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
m.contourf(Xm, Ym, V[-400/dz,...]*100, np.arange(-50, 50+2,2), cmap=plt.cm.RdBu_r, extend='both')
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

# plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/fm_vel_horiz.pdf')
# plt.savefig('/home/rsoutelino/prod/congressos/roms_wksp2012/fm_vel_horiz.png', dpi=300)
###############################################################################################


## schematic figure:#############################################################################

pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(V[:,0,:])*100,np.arange(-20,20+3,3), cmap=plt.cm.RdBu, alpha=1)
plt.colorbar(shrink=0.9, extend='both')
plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
plt.text(-42.6, 5, '$cm.s^{-1}$')
plt.text(-42.6, -1250, '$Vel$')

delta  = delta/111.0
x0     = -43.58
z0     = z0max.min()
height = 5   
offset = 0.05
plt.plot((-delta / 2.+x0, delta / 2.+x0), (z0, z0), 'k', linewidth = 2)
plt.text(x0+0.05,z0-70, r'$\delta$', {'color' : 'k', 'fontsize' : 24})
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0-height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0-height), 'k', linewidth = 2)

z0     = 0
plt.plot((-delta / 2.+x0, delta / 2.+x0), (z0, z0), 'k', linewidth = 2)
plt.text(x0+0.05,z0-70, r'$\delta$', {'color' : 'k', 'fontsize' : 24})
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0-height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0-height), 'k', linewidth = 2)

delta   = 2*np.sqrt(2*np.log1p(2)) * ds_NBUC
x0      = z0max.min()
z0      = -43.58
height  = 0.008
offset  = 20
plt.plot((z0, z0), (x0, delta / 2.+x0), '0.5', linewidth = 2)
plt.text(z0-0.16,x0+50, r'$\frac{1}{2}$', {'color' : '0.5', 'fontsize' : 13})
plt.text(z0-0.12,x0+50, r'$\delta_t$', {'color' : '0.5', 'fontsize' : 18})
plt.plot((z0, z0+height), (delta/2.+x0, delta/2.+x0-offset), '0.5', linewidth = 2)
plt.plot((z0, z0-height), (delta/2.+x0, delta/2.+x0-offset), '0.5', linewidth = 2)

delta   = 2*np.sqrt(2*np.log1p(2)) * db_NBUC
x0      = z0max.min()
z0      = -43.58
height  = 0.008
offset  = 20
plt.plot((z0, z0), (x0, -delta / 2.+x0), 'b', linewidth = 2)
plt.text(z0-0.16,x0-280, r'$\frac{1}{2}$', {'color' : 'b', 'fontsize' : 13})
plt.text(z0-0.12,x0-280, r'$\delta_b$', {'color' : 'b', 'fontsize' : 18})
plt.plot((z0, z0+height), (-delta/2.+x0, -delta/2.+x0+offset), 'b', linewidth = 2)
plt.plot((z0, z0-height), (-delta/2.+x0, -delta/2.+x0+offset), 'b', linewidth = 2)

x0     = -43.58
z0     = z0max.min()
plt.plot(x0,z0, 'ko',markersize=8)
plt.plot(x0,z0, 'wo',markersize=6)

delta   = 2*np.sqrt(2*np.log1p(2)) * d_BC
x0      = 0
z0      = -43.58
height  = 0.008
offset  = 20
plt.plot((z0, z0), (x0, -delta / 2.+x0), 'r', linewidth = 2)
plt.text(z0-0.16,x0-110, r'$\frac{1}{2}$', {'color' : 'r', 'fontsize' : 13})
plt.text(z0-0.12,x0-110, r'$\delta_b$', {'color' : 'r', 'fontsize' : 18})
plt.plot((z0, z0+height), (-delta/2.+x0, -delta/2.+x0+offset), 'r', linewidth = 2)
plt.plot((z0, z0-height), (-delta/2.+x0, -delta/2.+x0+offset), 'r', linewidth = 2)

x0     = -43.58
z0     = 0
plt.plot(x0,z0, 'ko',markersize=8)
plt.plot(x0,z0, 'wo',markersize=6)

plt.plot(x0,-1100, 'wv',markersize=6)
plt.text(x0-0.1,-1150, r'$x = x_c$', {'color' : 'k', 'fontsize' : 16})

plt.plot(-43,-500, 'w>',markersize=6)
plt.text(-43+0.07,-514, r'$z = z_c$', {'color' : 'k', 'fontsize' : 16})

plt.plot(-43,0, 'w>',markersize=6)
plt.text(-43+0.07,-14, r'$z = 0$', {'color' : 'k', 'fontsize' : 16})

plt.plot((-44.4, -42.7),(-1200,-1200), 'k--', linewidth=5, alpha=0.1)

plt.plot(-43.2,-1200, 'w>',markersize=6)
plt.text(-43.2+0.07,-1214, r'$z = 1200 m$', {'color' : 'k', 'fontsize' : 16})
plt.title('South')

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = (-44.3,-1150)
plt.text(ax, ay, "a", ha="center", va="center", size=12, bbox=bbox_props)

# plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/bc-nbuc_diagram_south.pdf')
# plt.savefig('/home/rsoutelino/prod/congressos/roms_wksp2012/bc-nbuc_diagram_south.png', dpi=300)
#==========================================================================

pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(V[:,-1,:])*100,np.arange(-50,50+6,6), cmap=plt.cm.RdBu, alpha=1)
plt.colorbar(shrink=0.9, extend='both')
plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
plt.text(-42.6, 5, '$cm.s^{-1}$')
plt.text(-42.6, -1250, '$Vel$')

delta  = 100/111.0
x0     = -43.58
z0     = z0max.max()
height = 5   
offset = 0.05
plt.plot((-delta / 2.+x0, delta / 2.+x0), (z0, z0), 'k', linewidth = 2)
plt.text(x0+0.05,z0-70, r'$\delta$', {'color' : 'k', 'fontsize' : 24})
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((-delta/2.+x0, -delta/2.+x0+offset), (z0, z0-height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0+height), 'k', linewidth = 2)
plt.plot((delta/2.+x0, delta/2.+x0-offset), (z0, z0-height), 'k', linewidth = 2)

delta   = 2*np.sqrt(2*np.log1p(2)) * ds_NBUC
x0      = z0max.max()
z0      = -43.58
height  = 0.008
offset  = 20
plt.plot((z0, z0), (x0, delta / 2.+x0), '0.5', linewidth = 2)
plt.text(z0-0.16,x0+50, r'$\frac{1}{2}$', {'color' : '0.5', 'fontsize' : 13})
plt.text(z0-0.12,x0+50, r'$\delta_t$', {'color' : '0.5', 'fontsize' : 18})
plt.plot((z0, z0+height), (delta/2.+x0, delta/2.+x0-offset), '0.5', linewidth = 2)
plt.plot((z0, z0-height), (delta/2.+x0, delta/2.+x0-offset), '0.5', linewidth = 2)

delta   = 2*np.sqrt(2*np.log1p(2)) * db_NBUC
x0      = z0max.max()
z0      = -43.58
height  = 0.008
offset  = 20
plt.plot((z0, z0), (x0, -delta / 2.+x0), 'b', linewidth = 2)
plt.text(z0-0.16,x0-280, r'$\frac{1}{2}$', {'color' : 'b', 'fontsize' : 13})
plt.text(z0-0.12,x0-280, r'$\delta_b$', {'color' : 'b', 'fontsize' : 18})
plt.plot((z0, z0+height), (-delta/2.+x0, -delta/2.+x0+offset), 'b', linewidth = 2)
plt.plot((z0, z0-height), (-delta/2.+x0, -delta/2.+x0+offset), 'b', linewidth = 2)

x0     = -43.58
z0     = z0max.max()
plt.plot(x0,z0, 'ko',markersize=8)
plt.plot(x0,z0, 'wo',markersize=6)


plt.plot(x0,-800, 'wv',markersize=6)
plt.text(x0-0.1,-850, r'$x = x_c$', {'color' : 'k', 'fontsize' : 16})

plt.plot(-43,-200, 'w>',markersize=6)
plt.text(-43+0.07,-214, r'$z = z_c$', {'color' : 'k', 'fontsize' : 16})

plt.plot(-43,0, 'w>',markersize=6)
plt.text(-43+0.07,-14, r'$z = 0$', {'color' : 'k', 'fontsize' : 16})

plt.plot((-44.4, -42.7),(-1200,-1200), 'k--', linewidth=5, alpha=0.1)

plt.plot(-43.2,-1200, 'w>',markersize=6)
plt.text(-43.2+0.07,-1214, r'$z = 1200 m$', {'color' : 'k', 'fontsize' : 16})

plt.title('North')

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = (-44.3,-1150)
plt.text(ax, ay, "b", ha="center", va="center", size=12, bbox=bbox_props)

# plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/bc-nbuc_diagram_north.pdf')
# plt.savefig('/home/rsoutelino/prod/congressos/roms_wksp2012/bc-nbuc_diagram_north.png', dpi=300)

#############################################  
   
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


###############################################################################################
pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(RHO[:,0,:])-1000,np.arange(23.2,27.8,0.2), cmap=plt.cm.Spectral_r, alpha=1)
plt.colorbar(shrink=0.9, extend='both')
plt.contour(X[0,:],Z, np.squeeze(RHO[:,0,:])-1000,np.arange(23.2,27.8,0.2), colors='0.4') 
plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
plt.plot(x0,-1100, 'wv',markersize=6)
plt.text(x0-0.1,-1150, r'$x = x_c$', {'color' : 'k', 'fontsize' : 16})
plt.plot(-43,-500, 'w>',markersize=6)
plt.text(-43+0.07,-514, r'$z = z_c$', {'color' : 'k', 'fontsize' : 16})
plt.plot(-43,0, 'w>',markersize=6)
plt.text(-43+0.07,-14, r'$z = 0$', {'color' : 'k', 'fontsize' : 16})
plt.plot((-44.4, -42.7),(-1200,-1200), 'k--', linewidth=5, alpha=0.1)
plt.plot(-43.2,-1200, 'w>',markersize=6)
plt.text(-43.2+0.07,-1214, r'$z = 1200 m$', {'color' : 'k', 'fontsize' : 16})
plt.text(-42.6, 5, '$kg.m^{-3}$')
plt.text(-42.6, -1250, '$\sigma_\Theta$')
plt.title('South')

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = (-44.3,-1150)
plt.text(ax, ay, "a", ha="center", va="center", size=12, bbox=bbox_props)

# plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/bc-nbuc_diagram_rho_south.pdf')
# plt.savefig('/home/rsoutelino/prod/congressos/roms_wksp2012/bc-nbuc_diagram_rho_south.png', dpi=300)

pp = plt.figure(figsize=(8,8), facecolor='w')
ax = plt.axes()
cc = plt.contourf(X[0,:],Z, np.squeeze(RHO[:,-1,:])-1000,np.arange(23.2,27.8,0.2), cmap=plt.cm.Spectral_r, alpha=1) 
plt.colorbar(shrink=0.9, extend='both')
plt.contour(X[0,:],Z, np.squeeze(RHO[:,-1,:])-1000,np.arange(23.2,27.8,0.2), colors='0.4') 
plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
plt.plot(x0,-800, 'wv',markersize=6)
plt.text(x0-0.1,-850, r'$x = x_c$', {'color' : 'k', 'fontsize' : 16})
plt.plot(-43,-200, 'w>',markersize=6)
plt.text(-43+0.07,-214, r'$z = z_c$', {'color' : 'k', 'fontsize' : 16})
plt.plot(-43,0, 'w>',markersize=6)
plt.text(-43+0.07,-14, r'$z = 0$', {'color' : 'k', 'fontsize' : 16})
plt.plot((-44.4, -42.7),(-1200,-1200), 'k--', linewidth=5, alpha=0.1)
plt.plot(-43.2,-1200, 'w>',markersize=6)
plt.text(-43.2+0.07,-1214, r'$z = 1200 m$', {'color' : 'k', 'fontsize' : 16})
plt.text(-42.6, 5, '$kg.m^{-3}$')
plt.text(-42.6, -1250, '$\sigma_\Theta$')
plt.title('North')

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = (-44.3,-1150)
plt.text(ax, ay, "b", ha="center", va="center", size=12, bbox=bbox_props)

# plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/bc-nbuc_diagram_rho_north.pdf')
# plt.savefig('/home/rsoutelino/prod/congressos/roms_wksp2012/bc-nbuc_diagram_rho_north.png', dpi=300)
stop

for k in range(22):
    plt.figure(figsize=(8,8), facecolor='w')
    plt.contourf(X[0,:],Z, np.squeeze(V[:,k,:])*100,np.arange(-30,30+3,3), cmap=plt.cm.RdBu, extend='both')    
    plt.colorbar(shrink=0.9, extend='both')
    plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
    plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
    plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
    plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
    plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
    plt.text(-42.6, 5, '$cm.s^{-1}$')
    plt.text(-42.6, -1250, '$Vel$')
    plt.savefig('frame%02i.png' %k)

for k in range(22):
    plt.figure(figsize=(8,8), facecolor='w')
    plt.contourf(X[0,:],Z, np.squeeze(RHO[:,k,:])-1000,np.arange(23.2,27.8,0.2), 
                 cmap=plt.cm.Spectral_r, extend='both') 
    plt.colorbar(shrink=0.9, extend='both')
    plt.contour(X[0,:],Z, np.squeeze(RHO[:,k,:])-1000,np.arange(23.2,27.8,0.2), colors='k') 
    plt.axis([-44.4, -42.7, -1250, 50]); plt.axis('off')
    plt.plot([-44.4, -44.4, -42.7, -42.7, -44.4],[-1250, 50, 50, -1250, -1250],'k', linewidth=2)
    plt.plot((-44.4, -42.7),(0,0), 'k--', linewidth=5, alpha=0.1)
    plt.text(-44.3, -1300, '$X$: Cross-stream axis $\longrightarrow$');
    plt.text(-44.5, -860, '$Z$: Vertical axis  $\longrightarrow$', rotation=90);
    plt.text(-42.6, 5, '$cm.s^{-1}$')
    plt.text(-42.6, -1250, '$Vel$')
    plt.savefig('rho%s.png' %k)







stop2
###############################################################################################
Xm, Ym = m(X, Y)  
Xm = Xm.transpose()
Ym = Ym.transpose()

for k in range(21):
  plt.figure(facecolor='w', figsize=(7,8))
  m.plot(Xm, Ym, 'k', linewidth=4, alpha=0.5)
  m.plot(Xm[:,k], Ym[:,k], 'r', linewidth=5 )
  m.contourf(lonm, latm, zb, colors=('0.9'), alpha=0.5)  
  m.drawcoastlines(zorder=5)
  m.drawcountries(zorder=4)
  m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
  m.drawparallels(np.arange(latr.min(), latr.max(), 2),
     labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
  m.drawmeridians(np.arange(lonr.min(), lonr.max(), 2),
     labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
  plt.savefig('transect%s.png' %k)

