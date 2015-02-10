#!/usr/bin/env python
# -*- coding:utf-8 -*-
###################################################################
# SCRIPT FOR VIZUALIZATION OF AVISO ALTIMETRY 
#      PRODUCT AND COMPUTING MEAN FIELDS
# - HORIZONTAL MAPS
# - Rafael Soutelino - rsoutelino@gmail.com
# - last update: August, 2010
#
###################################################################
import pylab as pl
#import numpy as np
#import matplotlib.pyplot as plt
import scipy.io as sp
from mpl_toolkits.basemap import Basemap
import datetime as dt
import numpy as np

m = Basemap(projection='cyl',llcrnrlat=-25,urcrnrlat=-9,
	llcrnrlon=-46,urcrnrlon=-30,lat_ts=0,resolution='l')

####################################################################
print " "
print "Loading AVISO data from matlab computations"
dat = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/aviso_omar.mat')
mag = dat.pop('mag')
lon = dat.pop('lon')
lat = dat.pop('lat')
um = dat.pop('um')
vm = dat.pop('vm')
tt = dat.pop('tt')
Vab = dat.pop('Vab')
Vcf = dat.pop('Vcf')

lon,lat = m(lon,lat)
 
# masking 
f = pl.find(pl.isnan(mag) == 1); 
a, b = mag.shape; 
mag = pl.ravel(mag); um = pl.ravel(um); vm = pl.ravel(vm);
mag[f] = 0; um[f] = 0; vm[f] = 0;
mag.shape = (a,b); um.shape = (a,b); vm.shape = (a,b)
mag = pl.ma.masked_where(mag == 0, mag)
um = pl.ma.masked_where(um == 0, um)
vm = pl.ma.masked_where(vm == 0, vm)

print " "
print "Loading etopo data ....... " 
etopo = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/etopo2_atlsul.mat')
xb = etopo.pop('xb')
yb = etopo.pop('yb')
zb = etopo.pop('zb')
xb,yb = m(xb,yb)
zbm = pl.ma.masked_where(zb <= -1000, zb)

print " "
print "Plotting figure ............." 
#####################################################################
fig1 = pl.figure(1,figsize=(16,8.5),facecolor='w')

print " "
print "Plotting Map ..............."
#####################################################################
p1 = pl.subplot(121)

pc = m.contourf(lon,lat,mag, np.arange(5, 21, 0.5), cmap=pl.cm.Reds, extend='both')
m.quiver(lon,lat,um,vm,scale=300,color='k')
m.contourf(xb,yb,zbm,colors=('w'),linewidth=0)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(zorder=3)
m.drawparallels(pl.arange(-24,0,3),
	labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(pl.arange(-45,-20,3),
	labels=[0, 0, 0, 1],dashes=[1,1000],zorder=7)
ax,ay = m(-37,-18); m.plot(ax,ay,'*y',mec='k',mfc='y',ms=15)
ax,ay = m(-40.3,-23); m.plot(ax,ay,'*y',mec='k',mfc='y',ms=15)
ax,ay = m(-41,-14.5); pl.text(ax,ay,u"Ilhéus",color='k',fontsize=12)
ax,ay = m(-40.3,-12.5); pl.text(ax,ay,"Salvador",color='k',fontsize=12)
ax,ay = m(-41.75,-17.6); pl.text(ax,ay,"Caravelas",color='k',fontsize=12)
ax,ay = m(-42.6,-21); pl.text(ax,ay,"Vitoria",color='k',fontsize=12)
ax,ay = m(-44,-11); pl.text(ax,ay,"BRASIL",color='k',
	fontsize=14,fontweight='bold')

ax = fig1.add_axes([0.06, 0.27, 0.01, 0.46])
cbar = pl.colorbar(pc,cax=ax,orientation='vertical',alpha=1)
cbar.set_label("cm s$^{-1}$")
p1.set_title('AVISO Geostrophic Velocities')

#########################################################################
print " " 
print "Plotting Abrolhos series .............."
p2 = pl.subplot(222)

xmin = int(tt[0]); xmax = int(tt[-1])

pl.plot([tt[0],tt[-1]],[0, 0],'k'); pl.grid()
pl.plot(tt,Vab,'k',linewidth=1)
pl.plot(tt,pl.mean(Vab)*pl.ones(pl.shape(tt)),'r',linewidth=2)
p2.set_xlim(xmin, xmax); p2.set_ylim(-50, 50)
p2.set_ylabel("V [cm s$^{-1}$]")
p2.set_title(r"Meridional Velocity off Abrolhos Bank - 18$^\circ$S")
p2.set_axis_bgcolor('0.95')
p2.xaxis.set_major_formatter(pl.DateFormatter('%b'))
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
p2.text(dt.datetime(2008,10,1),-40,r"Southward = 63%",
	ha="center", va="center", size=12,bbox=bbox_props)
p2.text(dt.datetime(2008,10,1),40,r"Mean Velocity = -1.8 cm s$^{-1}$",
	ha="center", va="center", size=12,bbox=bbox_props,color='r')
#

#########################################################################
print " "
print "Plotting Cabo Frio Series ............"
p3 = pl.subplot(224)

pl.plot([tt[0],tt[-1]],[0, 0],'k'); pl.grid()
pl.plot(tt,Vcf,'k',linewidth=1)
pl.plot(tt,pl.mean(Vcf)*pl.ones(pl.shape(tt)),'r',linewidth=2)
p3.set_xlim(xmin, xmax); p3.set_ylim(-50, 50)
p3.set_xlabel("Jan, 2007 - Aug, 2009")
p3.set_ylabel("V [cm s$^{-1}$]")
p3.set_title(r"Meridional Velocity off Cabo Frio - 22$^\circ$S")
p3.set_axis_bgcolor('0.95')
p3.xaxis.set_major_formatter(pl.DateFormatter('%b'))
p3.text(dt.datetime(2008,10,1),-40,r"Southward = 89%",
	ha="center", va="center", size=12,bbox=bbox_props)
p3.text(dt.datetime(2008,10,1),40,r"Mean Velocity = -10 cm s$^{-1}$",
	ha="center", va="center", size=12,bbox=bbox_props,color='r')


fig1.tight_layout(pad=3)
print " "
print "Writting PDF file ............"
pl.savefig('figures/aviso_time_series.pdf')
pl.show()


