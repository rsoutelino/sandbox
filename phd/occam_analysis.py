#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script for vizualization of OCCAM fields
#
# Rafael Soutelino - rsoutelino@gmail.com
#
#
# Last modification: Apr, 2011
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
from datetime import datetime

def near(x,x0):
	"""
	Find the index where x has the closer value to x0
	"""
	
	dx = x - x0
	dx = np.abs(dx)
	fn = np.where( dx == dx.min() )
	fn = fn[0][0]
	
	return fn

#####################################################################

m = Basemap(projection='merc',llcrnrlat=-23.0,urcrnrlat=-5.0,
	llcrnrlon=-41.0,urcrnrlon=-20.0,lat_ts=0,resolution='l')

# defining sub-area
dataset = nc.Dataset('data/jan2003.h5m1subvol')
lon = dataset.variables['LONGITUDE_T'][:]; lon = lon-360
lat = dataset.variables['LATITUDE_T'][:]
fb1 = near(lon, -33); fb2 = near(lon, -20); 
fa1 = near(lat, -20); fa2 = near(lat, -6);  

yeararray  = range(1988, 2005)
montharray = ('jan', 'feb', 'mar', 'apr', 'may', 'jun', 
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

# creating empty arrays
c=0; umean=0; vmean=0
secmax=[]
year=[]
month=[]
secmaxlon=[]; secmaxlat=[]

for Y in yeararray:
    c = 0
    for M in montharray:
        c = c + 1
        month = np.hstack((month, c))
        year = np.hstack((year, Y))
        dataset = nc.Dataset('data/'+ M + str(Y) +'.h5m1subvol')
        #lon = dataset.variables['LONGITUDE_T'][fb1:fb2]; lon = lon-360
        #lat = dataset.variables['LATITUDE_T'][fa1:fa2]; 
        lon = dataset.variables['LONGITUDE_T'][:]; lon = lon-360
        lat = dataset.variables['LATITUDE_T'][:];      
        
        #temp  = dataset.variables['POTENTIAL_TEMPERATURE__MEAN_'][:]
        #U     = dataset.variables['U_VELOCITY__MEAN_'][0,fa1:fa2,fb1:fb2]; U = U/100
        U     = dataset.variables['U_VELOCITY__MEAN_'][:]; U = U/100
        V     = dataset.variables['V_VELOCITY__MEAN_'][:]; V = V/100
        #depth = dataset.variables['DEPTH'][:]; depth = depth/100

        #etopo = sp.loadmat('/home/rsoutelino/rsoutelino/mestrado/proc/common/etopo2_leste.mat')
        #xb = etopo.pop('xb')
        #yb = etopo.pop('yb')
        #zb = etopo.pop('zb')
        #xb,yb = m(xb,yb)

        # masking topography to plot continental shelf mask
        #zb = np.ma.masked_where(zb <= -1000, zb)

        #d  = 2
        #p  = 60
        #pn = near(depth, p); 

        #u    = np.squeeze(U[pn,:,:]); 
        #v    = np.squeeze(V[pn,:,:]); 
        #temp = np.squeeze(temp[pn,:,:]); 
        #lont, latt = np.meshgrid(lont, latt)
        #lonu, latu = np.meshgrid(lonu, latu)
        #lontp,lattp = m(lont,latt)
        #lonup,latup = m(lonu,latu)

        #umean = umean + u
        #vmean = vmean + v
        lon, lat = np.meshgrid(lon, lat)
        
        #secmax = np.hstack((secmax, U.min()))
        #fmax = np.where(U == U.min())
        #secmaxlon = np.hstack((secmaxlon, lon[fmax]))
        #secmaxlat = np.hstack((secmaxlat, lat[fmax]))
        
        plt.figure()
        plt.quiver(lon[::5,::5], lat[::5,::5], np.squeeze(U[10,::5,::5]), np.squeeze(V[10,::5,::5]))
        plt.axis('equal')
        plt.axis([-42, -19, -24, -4])
        plt.savefig('lixo' + str(Y) + '_' + str(c) + '.png')
        


plt.figure()
plt.plot(secmax*100,'g'); plt.grid()
plt.show()

secmaxlonm, secmaxlatm = m(secmaxlon, secmaxlat) 

plt.figure()
m.plot(secmaxlonm, secmaxlatm, 'b.')
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-6,4),labels=[0, 1, 0, 0],dashes=[1,3],zorder=6)
m.drawmeridians(np.arange(-41,-20,4),labels=[0, 0, 0, 1],dashes=[1,3],zorder=7)
plt.show()































