#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  roms_analysis1.py -- Analysis of ROMS fields
#  - time-averaged horizontal fields and vertical sections
#
#  Last Modification: May, 2011
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


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
days        = np.arange(29,360) # length of the series

nk          = 400           # vertical level [m]
lonplot     = (-41, -33.7)  # longitude limits for plotting (OEII grid)
latplot     = (-22.5, -12)  # latitude limits for plotting
#lonplot     = (-41, -32)  # longitude limits for plotting (FULL GRID)
#latplot     = (-23, -10)  # latitude limits for plotting
hsc         = (26, 30)    # color scale limits for horizontal plots
d           = 4           # subsampling factor for decent quiver plot [integer]
sc          = 3           # scale to quiver arrows

vsc         = (-0.3, 0.3+0.03)
vst         = 0.03        # contour interval

niN         = -12         # latitude for NBUC section
niC         = -17         # latitude for RC Eddy
niS         = -19         # longitude for Abrolhos Eddy

XlimN       = (-38, -35)  # x axis limits for northern section
XlimC       = (-39, -35)  # x axis limits for southern section
XlimS       = (-39, -35)  # x axis limits for eastern section

Zlim        = (-1500, 0)  # z axis limits

# READING ROMS FIELDS ============================================================

# TO REMOVE LOOPING OVER DAYS:
#    comment next two lines, 
#    unident the whole loop
#    uncomment print statements
#    uncomment for statement to load variables
#    ident whatever pertains to the for statement
#    uncomment line that computes the averages
#    uncomment plt.show()
#    comment close('all')

#for days in np.arange(1,2,10):
    #print 'Day: '+ str(days); 

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
T=0; U=0; V=0; W=0;
for l in days: 
    T    = T + outfile.variables['temp'][l,...]
    U    = U + outfile.variables['u'][l,...]
    V    = V + outfile.variables['v'][l,...]
    #W    = W + outfile.variables['w'][l,:-1,...]

T = T / len(days); U = U / len(days); V = V / len(days); 
W = W / len(days);

km, im, jm = T.shape

print ' \n' + '==> ' + '  GETTING DEPTHS OF S-LEVELS ...\n' + ' '
zt   = get_depths(outfile,grdfile,0,'temp')
zu   = get_depths(outfile,grdfile,0,'u')
zv   = get_depths(outfile,grdfile,0,'v')

# horizontal ==========================================

print ' \n' + '==> ' + '  INTERPOLATING FROM S --> Z ...\n' + ' '
u = 0*lonu; v = 0*lonv; t = 0*lon; s = 0*lon; w = 0*lon

for a in range (0, im):
    for b in range(0, jm):
        t[a,b] = np.interp(-nk, zt[:, a, b], T[:, a, b] )
        #s[a,b] = np.interp(-nk, zt[:, a, b], S[:, a, b] )
        #w[a,b] = np.interp(-nk, zt[:, a, b], W[:, a, b] )

for a in range (0, im):
    for b in range(0, jm-1):
        u[a,b] = np.interp(-nk, zu[:, a, b], U[:, a, b] )

for a in range (0, im-1):
    for b in range(0, jm):
        v[a,b] = np.interp(-nk, zv[:, a, b], V[:, a, b] )
        
print ' \n' + '==> ' +\
        '  INTERPOLATING VECTORIAL VARIABLES TO SCALAR RHO POINTS ...\n' + ' '  

u   = griddata(lonu.ravel(), latu.ravel(), u.ravel(), lon, lat) 
v   = griddata(lonv.ravel(), latv.ravel(), v.ravel(), lon, lat)
eke = (1./2.) * (u**2 + v**2)

# masking in land values
t = np.ma.masked_where(h < nk, t)
#s = np.ma.masked_where(h < nk, s)
u   = np.ma.masked_where(h < nk, u)
v   = np.ma.masked_where(h < nk, v)
eke = np.ma.masked_where(h < nk, eke)
#w = np.ma.masked_where(h < nk, w)

print ' \n' + '==> ' + '  PLOTTING FIGURE ...\n' + ' '   

# setting up map projection
m = Basemap(projection='merc',
    llcrnrlat = latplot[0], urcrnrlat = latplot[1], 
    llcrnrlon = lonplot[0], urcrnrlon = lonplot[1],
      lat_ts=0, resolution='i')
      
lonm, latm = m(lon, lat)

up   = u[::d, ::d]
vp   = v[::d, ::d]
lonp = lonm[::d, ::d]
latp = latm[::d, ::d]
up   = np.ma.masked_where(up > 1e30, up)
vp   = np.ma.masked_where(vp > 1e30, vp)
eke  = np.ma.masked_where(eke > 1e30, eke)
#wp   = np.ma.masked_where(w > 1e30, w)

vmax = np.sqrt(up**2 + vp**2); vmax = int(vmax.max()*100)
days = days + 1

# vel over temp or vel over salt or vel over h -----------------------------
fig1 = plt.figure(1, figsize=(6,8), facecolor='w')

#m.pcolormesh(lonm,latm,prop, vmin=hsc[0], vmax=hsc[1],
    #cmap=plt.cm.Spectral_r)
#m.pcolormesh(lonm,latm,wp, vmin=-1e-4, vmax=1e-4,cmap=plt.cm.RdBu_r)
con = m.contourf(lonm,latm,h, 30,cmap=plt.cm.Blues,alpha=0.2)
#m.contourf(lonm,latm,t, 30,cmap=plt.cm.Spectral_r,alpha=0.8)
#m.contourf(lonm,latm, np.log10(eke*10000), 30, cmap=plt.cm.Spectral_r,alpha=0.3) 
#plt.colorbar()
#zb = np.ma.masked_where(h >= 1000, h)
#m.contourf(lonm, latm, zb, colors=('0.7'), alpha=0.5) 
m.quiver(lonp, latp, up*sc, vp*sc, scale=10)
text = str(vmax) +" cm/s";
ax,ay = m( lonplot[0] + 0.5, latplot[1] - 1.6 ) 
m.quiver(ax,ay,vmax/100.0*sc,0,scale=10,color='k',zorder=10, headwidth=1, headlength=2)
ax,ay = m(lonplot[0] + 0.5, latplot[1] - 1.5 ); 
plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latplot[0], latplot[1], 1),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonplot[0], lonplot[1], 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)

tx, ty = m( lonplot[0] + 0.5, latplot[1] - 0.5 )
if days.size > 1:
    text = 'Days: ' + str(days[0]) +' - '+ str(days[-1])
else: 
    text = 'Day: ' + str(days[0]) 

plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')
tx, ty = m( lonplot[0] + 0.5, latplot[1] - 1 )
text = 'Z: '+ str(nk) +' m' 
plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')

bbox_props = dict(boxstyle="round", fc="w", ec="w", alpha=0.7)
tx, ty = m(-39,-18.5); plt.text(tx,ty,'AB')
tx, ty = m(-38.8,-16.2); plt.text(tx,ty,'RCB')
tx, ty = m(-36,-21); plt.text(tx,ty,'VTR')
#tx, ty = m(-39.2,-22); plt.text(tx,ty, 'BC', ha="center", va="center", bbox=bbox_props)
tx, ty = m(-37,-12.6); plt.text(tx,ty, 'NBUC', ha="center", va="center", bbox=bbox_props)
tx, ty = m(-39,-22); plt.text(tx,ty, 'NBUC', ha="center", va="center", bbox=bbox_props)
#tx, ty = m(-37,-17); plt.text(tx,ty, 'NBUC', ha="center", va="center", bbox=bbox_props)
#tx, ty = m(-36.5,-14.5); plt.text(tx,ty, 'BiSEC', ha="center", va="center", bbox=bbox_props)
#tx, ty = m(-37.5,-15); plt.text(tx,ty, 'IE', ha="center", va="center", bbox=bbox_props)
#tx, ty = m(-36.5,-17); plt.text(tx,ty, 'RCE', ha="center", va="center", bbox=bbox_props)
#tx, ty = m(-36.5,-19.2); plt.text(tx,ty, 'AE', ha="center", va="center", bbox=bbox_props)

ax = fig1.add_axes([0.2, 0.45, 0.02, 0.25])
cbar = plt.colorbar(con,cax=ax,orientation='vertical',alpha=0.3)
cbar.set_label('$Topog$ [m]')

plt.show()

if days.size > 1:
    plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/'+ expt +'_vel_day'+ 
      str(days[0]) +'-'+ str(days[-1]) +'_'+ str(nk) +'m.pdf')
else:
    plt.savefig('/home/rsoutelino/rsoutelino/prod/csr_phd/figures/'+ expt +'_vel_day'+ 
      str(days[0]) +'_'+ str(nk) +'m.pdf')


#plt.close('all')




STOP

# vel over EKE --------------------------------------------------------
fig5 = plt.figure(5, figsize=(6,8), facecolor='w')

m.pcolormesh(lonm,latm, np.log10(eke*10000), vmin=-3, vmax=3, cmap=plt.cm.jet); 
plt.colorbar(shrink=0.8)
#m.quiver(lonp, latp, up*sc, vp*sc, scale=10)
zb = np.ma.masked_where(h >= 100, h)
m.contourf(lonm, latm, zb, colors=('0.7'), alpha=0.5) 
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(latplot[0], latplot[1], 1),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(lonplot[0], lonplot[1], 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
tx, ty = m( lonplot[0] + 0.5, latplot[1] - 1 )
text = 'Day: ' + str(days[0]) +' - '+ str(days[-1])
plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')
tx, ty = m( lonplot[0] + 0.5, latplot[1] - 2 )
text = 'Z: '+ str(nk) +' m' 
plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')
plt.title('$log(EKe)$ [$cm^2 s^{-2}$]')

plt.show()
plt.savefig('figures/'+ expt +'/eke_day'+ 
   str(days[0]) +'-'+ str(days[-1]) +'_'+ str(nk) +'m.png')

# vertical ==============================================

f    = latv[:,0]    
fsec = near(f,niN)
xN   = lonv[0,:]
zN   = np.squeeze( zv[:, fsec, :] )
vN   = np.squeeze( V[: ,fsec , :] ) 

fsec = near(f,niS)
xS   = lonv[0,:]
zS   = np.squeeze( zv[:, fsec, :] )
vS   = np.squeeze( V[: ,fsec , :] ) 

fsec = near(f,niC)
xC   = lonv[0,:]
zC   = np.squeeze( zv[:, fsec, :] )
vC   = np.squeeze( V[: ,fsec , :] ) 

xS.shape = (1, xS.size)
xN.shape = (1, xN.size)
xC.shape = (1, xC.size)

xS = np.repeat(xS, roms.klevels, axis=0)
xN = np.repeat(xN, roms.klevels, axis=0)
xC = np.repeat(xC, roms.klevels, axis=0)

vS = np.ma.masked_where(vS > 1e10, vS)
vN = np.ma.masked_where(vN > 1e10, vN)
vC = np.ma.masked_where(vC > 1e10, vC)

fig2 = plt.figure(2, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xN,zN,vN, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimN[0], XlimN[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niN)) +'$^\circ$S')

fig3 = plt.figure(3, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xC,zC,vC, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimC[0], XlimC[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niC)) +'$^\circ$S')

fig4 = plt.figure(4, figsize=(12,6), facecolor='w')
ax = plt.axes(axisbg=('0.7'))
plt.contourf(xS,zS,vS, np.arange(vsc[0], vsc[1], vst),cmap=plt.cm.RdBu, extend='both')
plt.colorbar(); plt.axis([XlimS[0], XlimS[1], Zlim[0], Zlim[1]])
plt.xlabel('Longitude')
plt.ylabel('Depth')
plt.title('Days: ' + str(days[0]+1) +' - '+ str(days[-1]+1) +
   ' Time-averaged Cross-section Velocity @ '+ str(np.abs(niS)) +'$^\circ$S')

plt.show()

#fig4 = plt.figure(3, figsize=(12,6), facecolor='w')
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
    
    
