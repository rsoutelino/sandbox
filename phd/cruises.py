#!/usr/bin/env python
# -*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.io as sp
from mpl_toolkits.basemap import Basemap
import datetime as dt
import numpy as np
import matplotlib
import romslab

matplotlib.rcParams.update({'font.size': 10})
####################################################################
lims = [-41, -33.7, -22.5, -12]
m = Basemap(projection='cyl',llcrnrlat=-22.5,urcrnrlat=-12.0,
    llcrnrlon=-41.0,urcrnrlon=-33.7,lat_ts=0,resolution='l')

####################################################################
class Leste1:
    stations = np.loadtxt('/home/rsoutelino/mestrado/proc/hidrografia/posicoes_leste1.dat')
    lon, lat = stations[:,2], stations[:,1]
    lon = np.ma.masked_where(lat > -12.5, lon)
    lat = np.ma.masked_where(lat > -12.5, lat)

class Leste2:
    track = sp.loadmat('/home/rsoutelino/phd/data/leste2/adcp/leste2_xy.mat')
    xyt = track['xyt']
    lon, lat = xyt[0,:] - 360, xyt[1,:]
    lon = np.ma.masked_where(lat > -12.5, lon)
    lat = np.ma.masked_where(lat > -12.5, lat)

class ProAbro:
    track = sp.loadmat('/home/rsoutelino/phd/data/proab1/adcp/proab1_xy.mat')
    xyt = track['xyt']
    lon, lat = xyt[0,:] - 360, xyt[1,:] 
    lon = np.ma.masked_where(lat > -12.5, lon)
    lat = np.ma.masked_where(lat > -12.5, lat)
####################################################################

leste1, leste2, proab = Leste1(), Leste2(), ProAbro()

for dataset in [leste1, leste2, proab]:
    dataset.mlon, dataset.mlat = m(dataset.lon, dataset.lat)

print " "
print "Loading etopo data ....... " 
etopo = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/etopo2_atlsul.mat')
xb = etopo.pop('xb')
yb = etopo.pop('yb')
zb = etopo.pop('zb')
xb, yb, zb = romslab.subset(xb, yb, zb, lims[0], lims[1], lims[2], lims[3])
xb,yb = m(xb,yb)


fig1 = plt.figure(1,figsize=(15, 8),facecolor='w')

gs =gridspec.GridSpec(1, 3)

p1 = plt.subplot(gs[0, 0])
m.contourf(xb, yb, zb, np.arange(-5000, 0, 500), cmap=plt.cm.Blues_r)
m.plot(leste1.lon, leste1.lat, '*y', markersize=8)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(zorder=3)
m.drawparallels(np.arange(-24,0,2),
	labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-45,-20,2),
	labels=[0, 0, 0, 1],dashes=[1,1000],zorder=7)
p1.set_title('Oceano Leste I')

p2 = plt.subplot(gs[0, 1])
m.contourf(xb, yb, zb, np.arange(-5000, 0, 500), cmap=plt.cm.Blues_r)
m.plot(leste2.lon, leste2.lat, 'y', linewidth=2)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(zorder=3)
# m.drawparallels(np.arange(-24,0,2),
    # labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-45,-20,2),
    labels=[0, 0, 0, 1],dashes=[1,1000],zorder=7)
p2.set_title('Oceano Leste II')

p3 = plt.subplot(gs[0, 2])
m.contourf(xb, yb, zb, np.arange(-5000, 0, 500), cmap=plt.cm.Blues_r)
m.plot(proab.lon, proab.lat, 'y', linewidth=2)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(zorder=3)
# m.drawparallels(np.arange(-24,0,2),
    # labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-45,-20,2),
    labels=[0, 0, 0, 1],dashes=[1,1000],zorder=7)
p3.set_title('PRO-ABROLHOS')

fig1.tight_layout(pad=3)
plt.show()

print " "
print "Writting PDF file ............"
plt.savefig('figures/cruises.pdf')
  


