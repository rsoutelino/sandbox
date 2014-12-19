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

#####################################################################

m = Basemap(projection='merc',llcrnrlat=-23.0,urcrnrlat=-10.0,
	llcrnrlon=-41.0,urcrnrlon=-32.0,lat_ts=0,resolution='l')

y = 2003
mon = 'feb'


c = 0; umean = 0; vmean = 0
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

u = umean; v = vmean
#################################################
s = 2; # quiver scale

plt.figure(figsize=(6,8),facecolor='w')
# m.pcolormesh(lont[::d,::d], latt[::d,::d],temp[::d,::d])
m.quiver(lonup[::d,::d], latup[::d,::d], u[::d,::d]*s,v[::d,::d]*s,scale=10);
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
plt.title('OCCAM '+ mon + str(y) +' - Velocities') 

plt.show()

# plt.savefig('figures/occam_'+ mon + str(y) +'.png')


matdict = {'lonSoccam':lonSoccam, 'latSoccam':latSoccam, 'lonNoccam':lonNoccam, 'latNoccam':latNoccam, 'lonEoccam':lonEoccam, 'latEoccam':latEoccam, 'tempSoccam':tempSoccam, 'tempNoccam':tempNoccam, 'tempEoccam':tempEoccam, 'saltSoccam':saltSoccam, 'saltNoccam':saltNoccam, 'saltEoccam':saltEoccam}
