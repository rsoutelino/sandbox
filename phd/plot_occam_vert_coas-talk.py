#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script for vizualization of OCCAM fields
#
# Rafael Soutelino - rsoutelino@gmail.com
#
#
# Last modification: Oct, 2010
#####################################################################

print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 

# IMPORTING MODULES #################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
# import datetime as dt
import netCDF4 as nc
import scipy.io as sp

def near(x,x0):
	"""
	Find the index where x has the closer value to x0
	"""
	
	dx = x - x0
	dx = np.abs(dx)
	fn = np.where( dx == dx.min() )
	
	return fn

def transp1(lon, z, v):
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


def transp2(lat, z, u):
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

#####################################################################

m = Basemap(projection='merc',llcrnrlat=-23.0,urcrnrlat=-10.0,
	llcrnrlon=-41.0,urcrnrlon=-32.0,lat_ts=0,resolution='l')

#y = 2003
mon = 'y'


c = 0; umean = 0; vmean = 0; uu = 0; vv = 0
for y in range(2003, 2004):
    c = c + 1
    dataset = nc.Dataset('data/'+ mon + str(y) +'.h5m1subvol')
    lont = dataset.variables['LONGITUDE_T'][:]; lont = lont-360
    latt = dataset.variables['LATITUDE_T'][:]; 
    lonu = dataset.variables['LONGITUDE_U'][:]; lonu = lonu-360
    latu = dataset.variables['LATITUDE_U'][:]
    temp  = dataset.variables['POTENTIAL_TEMPERATURE__MEAN_'][:]
    U     = dataset.variables['U_VELOCITY__MEAN_'][:]; U = U/100
    V     = dataset.variables['V_VELOCITY__MEAN_'][:]; V = V/100
    depth = dataset.variables['DEPTH'][:]; depth = depth/100

    etopo = sp.loadmat('/home/rsoutelino/rsoutelino/mestrado/proc/common/etopo2_leste.mat')
    xb = etopo.pop('xb')
    yb = etopo.pop('yb')
    zb = etopo.pop('zb')
    xb,yb = m(xb,yb)

    # masking topography to plot continental shelf mask
    zb = np.ma.masked_where(zb <= -1000, zb)

    d  = 2
    p  = 60
    pn = near(depth, p); 

    u    = np.squeeze(U[pn,:,:]); 
    v    = np.squeeze(V[pn,:,:]); 
    temp = np.squeeze(temp[pn,:,:]); 
    lont, latt = np.meshgrid(lont, latt)
    lonu, latu = np.meshgrid(lonu, latu)
    lontp,lattp = m(lont,latt)
    lonup,latup = m(lonu,latu)

    umean = umean + u
    vmean = vmean + v
    uu = uu + U
    vv = vv + V
    

u = umean/c; v = vmean/c
U = uu/c; V = vv/c

#################################################
fS = 5; fN = -68; fE = -153
f1000 = 35

vst         = 5        # contour interval
vsc         = (-60, 60+vst) # color scale limits for vertical plot

s = 2; # quiver scale

plt.figure(1,facecolor='w')

#p1 = plt.subplot(131)
m.quiver(lonup[::d,::d], latup[::d,::d], u[::d,::d]*s,v[::d,::d]*s,scale=10);
m.plot(lonup[fN,:],latup[fN,:],'y',linewidth=10,alpha=0.5)
m.plot(lonup[fS,:],latup[fS,:],'g',linewidth=10,alpha=0.5)
#m.plot(lonup[:,fE],latup[:,fE],'y',linewidth=10,alpha=0.5)
m.contourf(xb,yb,zb,colors=('0.7'), alpha=0.5);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[1, 0, 0, 0],dashes=[1,3],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 0, 1],dashes=[1,3],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
#plt.title('OCCAM '+ mon + str(y) +' - Velocities') 

plt.figure(2,facecolor='w')

x,z = np.meshgrid(lonu[fN,:fE],-depth)
#TN = transp1(x,z,np.squeeze(V[:f1000,fN,:fE]))
p2 = plt.subplot(121,axis_bgcolor='y',alpha=0.5)
con2 = p2.pcolormesh(x,z, np.squeeze(V[:,fN,:fE]*100),vmin=-20, vmax=20, cmap=plt.cm.RdBu)
p2.contour(x,z, np.squeeze(V[:,fN,:fE]*100), np.arange(0, 100, 5), colors='w')
p2.contour(x,z, np.squeeze(V[:,fN,:fE]*100), np.arange(-100, 0, 1), colors='k')
p2.set_xlim(-36,-34); p2.set_ylim(-4000,0); 
p2.set_yticklabels([])
#p2.text(-34,-200,str(TN.round()) + ' Sv')
#cbar = plt.colorbar(con2,orientation='vertical',aspect=20)

#vst         = 1 # contour interval
#vsc         = (-10, 10+vst) # color scale limits for vertical plot
#y,z = np.meshgrid(latu[fS:fN,fE],-depth[:f1000])
#TE = transp2(y,z,np.squeeze(U[:f1000,fS:fN,fE]))
#p3 = plt.subplot(324,axis_bgcolor='y',alpha=0.5)
#con3 = p3.contourf(y,z, np.squeeze(U[:f1000,fS:fN,fE]*100),
    #np.arange(vsc[0], vsc[1], vst), cmap=plt.cm.RdBu)
#p3.set_xlim(-22.5,-10.7); p3.set_ylim(-1000,0);
#p3.text(-13,-800,str(TE.round()) + ' Sv')
#cbar = plt.colorbar(con3,orientation='vertical',aspect=20)

vst         = 2        # contour interval
vsc         = (-35, 35+vst) # color scale limits for vertical plot
x,z = np.meshgrid(lonu[fS,:fE],-depth)
#TS = transp1(x,z,np.squeeze(V[:,fS,:fE]))
p3 = plt.subplot(122,axis_bgcolor='g',alpha=0.5)
con3 = p3.pcolormesh(x,z, np.squeeze(V[:,fS,:fE]*100),vmin=-20, vmax=20, cmap=plt.cm.RdBu)
p3.contour(x,z, np.squeeze(V[:,fS,:fE]*100), np.arange(0, 100, 1), colors='w')
p3.contour(x,z, np.squeeze(V[:,fS,:fE]*100), np.arange(-100, 0, 1), colors='k')
p3.set_xlim(-41.5,-38); p3.set_ylim(-4000,0);
#p3.text(-35,-200,str(TS.round()) + ' Sv')
#cbar = plt.colorbar(con3,orientation='vertical',aspect=20, shrink=0.5)

plt.show()
###########################################################
stop

plt.figure(2,figsize=(14,8),facecolor='w')
p1 = plt.subplot(121)
m.quiver(lonup[::d,::d], latup[::d,::d], u[::d,::d]*s,v[::d,::d]*s,scale=10);
m.plot(lonup[fN,:],latup[fN,:],'g',linewidth=10,alpha=0.5)
m.plot(lonup[fS,:],latup[fS,:],'g',linewidth=10,alpha=0.5)
m.plot(lonup[:,fE],latup[:,fE],'y',linewidth=10,alpha=0.5)
m.contourf(xb,yb,zb,colors=('0.7'), alpha=0.5);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[0, 1, 0, 0],dashes=[1,3],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 0, 1],dashes=[1,3],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
plt.title('OCCAM '+ mon + '2003 - Velocities') 

x,z = np.meshgrid(lonu[fN,:fE],-depth[f1000:])
TN = transp1(x,z,np.squeeze(V[f1000:,fN,:fE]))
p2 = plt.subplot(322,axis_bgcolor='g',alpha=0.5)
con2 = p2.contourf(x,z, np.squeeze(V[f1000:,fN,:fE]*100),
    np.arange(vsc[0], vsc[1], vst), cmap=plt.cm.RdBu)
p2.set_xlim(-37,-33); p2.set_ylim(-5000,-1000);
p2.text(-34,-1500,str(TN.round()) + ' Sv')
cbar = plt.colorbar(con2,orientation='vertical',aspect=20)

vst         = 1 # contour interval
vsc         = (-10, 10+vst) # color scale limits for vertical plot
y,z = np.meshgrid(latu[fS:fN,fE],-depth[f1000:])
TE = transp2(y,z,np.squeeze(U[f1000:,fS:fN,fE]))
p3 = plt.subplot(324,axis_bgcolor='y',alpha=0.5)
con3 = p3.contourf(y,z, np.squeeze(U[f1000:,fS:fN,fE]*100),
    np.arange(vsc[0], vsc[1], vst), cmap=plt.cm.RdBu)
p3.set_xlim(-22.5,-10.7); p3.set_ylim(-5000,-1000);
p3.text(-13,-1500,str(TE.round()) + ' Sv')
cbar = plt.colorbar(con3,orientation='vertical',aspect=20)

vst         = 1        # contour interval
vsc         = (-10, 10+vst) # color scale limits for vertical plot
x,z = np.meshgrid(lonu[fS,:fE],-depth[f1000:])
TS = transp1(x,z,np.squeeze(V[f1000:,fS,:fE]))
p4 = plt.subplot(326,axis_bgcolor='g',alpha=0.5)
con4 = p4.contourf(x,z, np.squeeze(V[f1000:,fS,:fE]*100),
    np.arange(vsc[0], vsc[1], vst), cmap=plt.cm.RdBu)
p4.set_xlim(-42,-33); p4.set_ylim(-5000,-1000);
p4.text(-35,-1500,str(TS.round()) + ' Sv')
cbar = plt.colorbar(con4,orientation='vertical',aspect=20)


plt.show()
