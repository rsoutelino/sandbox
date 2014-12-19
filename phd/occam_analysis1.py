#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script for analysis of OCCAM fields
#
# FEATURES:
# - time-dependant SEC bifurcation based on monthly fields
# - time-dependant BC and NBUC transports
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
depth = dataset.variables['DEPTH'][:]; depth = depth/100
fj1 = near(lon, -41); fj2 = near(lon, -33); 
fi1 = near(lat, -20); fi2 = near(lat, -6);  
fk  = near(depth, 100)

yeararray  = range(1988, 2005)
montharray = ('jan', 'feb', 'mar', 'apr', 'may', 'jun', 
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

# creating empty arrays
c=0; umean=0; vmean=0
secmax=[]
years=[]
months=[]
seclat=[]
TBC=[]; TNBUC=[]

for Y in yeararray:
    c = 0
    for M in montharray:
        c = c + 1
        dataset = nc.Dataset('data/'+ M + str(Y) +'.h5m1subvol')
        lon = dataset.variables['LONGITUDE_T'][fj1:fj2]; lon = lon-360
        lat = dataset.variables['LATITUDE_T'][fi1:fi2];      
        U     = dataset.variables['U_VELOCITY__MEAN_'][:fk,fi1:fi2,fj1:fj2]; U = U/100
        V     = dataset.variables['V_VELOCITY__MEAN_'][:fk,fi1:fi2,fj1:fj2]; V = V/100
        depth = dataset.variables['DEPTH'][:fk]; depth = depth/100
        lon, lat = np.meshgrid(lon, lat)
        
        # SEC Position
        vsec = V.mean(axis=2)
        vsec = vsec.mean(axis=0)
        fsec = near(vsec, 0)
        seclat = np.hstack((seclat, lat[fsec,0]))
        months = np.hstack((months, c))
        years = np.hstack((years, Y))
        
        # BC Transport @ south of BiSEC
        Vbc = V[:,:10,:]; Vbc = Vbc.mean(axis=1)
        dx = np.diff(lon[0,:]); dx = dx.sum()*111000 # valid only for low latitudes
        dz = np.diff(depth); dz = dz.sum()
        Tbc = Vbc.mean() * dx * dz; Tbc = Tbc / 1000000
        TBC = np.hstack((TBC, Tbc))
        
        # NBUC Transport @ north of BiSEC
        Vnbuc = V[:,-10:,:]; Vnbuc = Vnbuc.mean(axis=1)
        Tnbuc = Vnbuc.mean() * dx * dz; Tnbuc = Tnbuc / 1000000 
        TNBUC = np.hstack((TNBUC, Tnbuc))
        


months = np.array(months, dtype=int)
years = np.array(years, dtype=int)
days  = months*0 + 15
time = []

for year, month, day in zip(years, months, days):
    time.append( datetime(year, month, day) )  


figure(1,facecolor='w', figsize=(16,8))

subplot(311)
plot(time, TNBUC,'b'); grid()
#xlabel('OCCAM Time')
ylabel('NBUC Transport')

subplot(312)
plot(time, seclat,'k'); grid()
#xlabel('OCCAM Time')
ylabel('BiSEC Latitude')

subplot(313)
plot(time, np.abs(TBC),'r'); grid()
xlabel('OCCAM Time')
ylabel('BC Transport')


























