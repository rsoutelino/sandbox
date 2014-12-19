#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script for vizualization of OCCAM fields
#
# Rafael Soutelino - rsoutelino@gmail.com
#
# FEATURES:
# - Mapping impinging agulhas rings and its effects on the western boundary
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
    print str(Y)
    for M in montharray:
        print str(M)
        c = c + 1
        month = np.hstack((month, c))
        year = np.hstack((year, Y))
        dataset1 = nc.Dataset('data/near-surface_big-domain/part1/'+ M + str(Y) +'.h5m1subvol')
        lon1  = dataset1.variables['LONGITUDE_T'][:]; lon1 = lon1-360
        lat1  = dataset1.variables['LATITUDE_T'][:];
        zeta1 = dataset1.variables['SEA_SURFACE_HEIGHT__MEAN_'][:]
        U1     = dataset1.variables['U_VELOCITY__MEAN_'][:]; U1 = U1/100
        V1     = dataset1.variables['V_VELOCITY__MEAN_'][:]; V1 = V1/100
        
        dataset2 = nc.Dataset('data/near-surface_big-domain/part2/'+ M + str(Y) +'.h5m1subvol')
        lon2  = dataset2.variables['LONGITUDE_T'][:]; 
        lat2  = dataset2.variables['LATITUDE_T'][:];
        zeta2 = dataset2.variables['SEA_SURFACE_HEIGHT__MEAN_'][:]
        #U2     = dataset2.variables['U_VELOCITY__MEAN_'][:]; U2 = U2/100
        #V2     = dataset2.variables['V_VELOCITY__MEAN_'][:]; V2 = V2/100       

        lon1, lat1 = np.meshgrid(lon1, lat1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        
        lon  = np.hstack( (lon1,  lon2) )
        lat  = np.hstack( (lat1,  lat2) )
        zeta = np.hstack( (zeta1, zeta2))
        
        plt.figure(facecolor='w', figsize=(14,8))
        
        m = Basemap(projection='merc',llcrnrlat=-40.0,urcrnrlat=-5.0,
            llcrnrlon=-50.0,urcrnrlon=25.0,lat_ts=0,resolution='l')
        lonm, latm = m(lon, lat)
        plt.axes([0.18, 0.56, 0.7, 0.4], axisbg='w')
        m.contourf(lonm, latm, np.squeeze(zeta), np.arange(-20, 21, 2), cmap=plt.cm.RdBu_r, extend='both'); # colorbar(); 
        m.drawcoastlines(zorder=5)
        #m.drawcountries(zorder=4)
        m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
        m.drawparallels((-25,-10),labels=[1, 1, 0, 0],dashes=[1,3],zorder=6)
        m.drawmeridians((-45,-25,5),labels=[0, 0, 1, 0],dashes=[1,3],zorder=7)
        x, y = m((-45, -25, -25, -45, -45),(-25, -25, -10, -10, -25))
        m.plot(x, y, 'k', zorder=10)
              
        m = Basemap(projection='merc',llcrnrlat=-25.0,urcrnrlat=-10.0,
            llcrnrlon=-45.0,urcrnrlon=-25.0,lat_ts=0,resolution='l') 
        lon1m, lat1m = m(lon1, lat1)
        plt.axes([0.47, 0.06, 0.4, 0.46], axisbg='w')
        m.contourf(lon1m, lat1m, np.squeeze(zeta1), np.arange(-5, 21, 1), cmap=plt.cm.RdBu_r, extend='both'); # colorbar(); 
        m.drawcoastlines(zorder=5)
        #m.drawcountries(zorder=4)
        m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
        m.drawparallels((-22, -14),labels=[0, 1, 0, 0],dashes=[1,3],zorder=6)
        m.drawmeridians((-41, -33),labels=[0, 0, 0, 1],dashes=[1,3],zorder=7)
        x, y = m((-41, -33, -33, -41, -41),(-21, -21, -13, -13, -21))
        m.plot(x, y, 'k', zorder=10)
               
        m = Basemap(projection='merc',llcrnrlat=-21.0,urcrnrlat=-13.0,
            llcrnrlon=-41.0,urcrnrlon=-33.0,lat_ts=0,resolution='l')
        lon1m, lat1m = m(lon1, lat1)
        plt.axes([0.2, 0.06, 0.3, 0.46], axisbg='w')
        m.contourf(lon1m, lat1m, np.squeeze(zeta1),  np.arange(-5, 21, 1), cmap=plt.cm.RdBu_r, alpha=0.5, extend='both')
        m.quiver(lon1m[::3,::3], lat1m[::3,::3], np.squeeze(U1[0,::3,::3]), np.squeeze(V1[0,::3,::3]), scale=5)
        m.drawcoastlines(zorder=5)
        #m.drawcountries(zorder=4)
        m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
        m.drawparallels(np.arange(-60,0,2),labels=[1, 0, 0, 0],dashes=[1,3],zorder=6)
        m.drawmeridians(np.arange(-60,30,2),labels=[0, 0, 0, 1],dashes=[1,3],zorder=7)
        
        plt.savefig('figures/agulhas_eddies' + str(Y) + '_' + str(c) + '.png')

        plt.close('all')
































