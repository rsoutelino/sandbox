#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  plot_roms.py -- Visualization of ROMS fields
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Jun, 2013
###################################################################################

import numpy as np
import sys
import scipy.io as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
# import netCDF4 as nc

plt.close('all')
# PYTHON SCRIPT START ######################################################

# Input parameters
expt = 'phd18'
zlev = 50
lims = [-40.0, -36, -19.5, -14]
d    = 3  # subsampling
s    = 5  # quiver size
sc   = 0.5 # quiver scale
spacer  = 0.5 # spacer between scale arrow and scale text 
 


if expt == 'phd15':
    panel = '(a)'
elif expt == 'phd18':
    panel = '(c)'
elif expt == 'phd16':
    panel = '(b)'


if expt == 'phd16':
    h = sp.loadmat('figures/h2.mat')['h']
else:
    snap = sp.loadmat('figures/day2555_50m.mat')
    h  = snap['h']

print ' \n' + '==> ' + '  READING FIELD ...\n' + ' '
avg  = sp.loadmat("averages/%s_run-average-vels_%sm.mat" %(expt, zlev))
lon  = avg['lon']
lat  = avg['lat']
u    = avg['umean']
v    = avg['vmean']

u = np.ma.masked_where(np.abs(u) > 2, u)
v = np.ma.masked_where(np.abs(v) > 2, v)


m = Basemap(projection='cyl', llcrnrlon = lims[0], urcrnrlon = lims[1], 
            llcrnrlat = lims[2], urcrnrlat = lims[3],
            lat_ts=0, resolution='i')


zb = np.ma.masked_where(h >= 100, h)
# u = np.ma.masked_where(h < 100, u)
# v = np.ma.masked_where(h < 100, v)


# plt.figure(figsize=(10,10), facecolor='w')
# m.quiver(lon[::d,::d], lat[::d,::d], u[::d,::d]*s, v[::d,::d]*s, scale=10)
# m.fillcontinents()
# m.drawcoastlines()
# m.drawmeridians(np.arange(lims[0], lims[1], 2),
#      labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
# m.drawparallels(np.arange(lims[2], lims[3], 1),
#      labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
# # scale
# ax, ay = m(lims[0]+spacer, lims[3]-spacer)
# m.quiver(ax, ay, sc, 0, scale=10, zorder=9)
# plt.text(ax, ay-spacer, "%s m/s" %(sc), color='k',fontsize=10,fontweight='bold')
# plt.title("Run-averaged velocities: %s m" %zlev)
# figname = "figures/run-average-vels_%s_%sm.png" %(expt, zlev)
# plt.savefig(figname)
# plt.show()


speed = np.sqrt(u*u + v*v)
# maxspeed = speed.max()
# speed = speed + 0.1
# speed[np.where(speed > 0.2)] = 0.2
d = 1
lw = 10*speed/speed.max()
print ' \n' + '==> ' + '  COMPUTING STREAMLINES ...\n' + ' '
fig = plt.figure(figsize=(8,10), facecolor='w')
splot = plt.streamplot(lon, lat, u, v, density=6, color='k', linewidth=lw, 
                                       arrowsize=3, arrowstyle='Fancy' )
# splot = plt.pcolormesh(lon, lat, speed, cmap=plt.cm.gray_r)
m.contourf(lon, lat, zb, colors=('0.9'), alpha=0.9, zorder=9)
m.fillcontinents(color='0.7', zorder=10)
m.drawcoastlines(zorder=11)
m.drawmeridians(np.arange(lims[0], lims[1], 1),
     labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
m.drawparallels(np.arange(lims[2], lims[3], 1),
     labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)

plt.text(-39.8, -17.5, panel, fontsize=16, fontweight='bold',zorder=12)
plt.text(-39.8, -14.5, "Max. Vel.:\n%0.2f m s$^{-1}$" %speed.max(), 
                       fontsize=12, fontweight='bold',zorder=12)

# ax2 = fig.add_axes([0.15, 0.44, 0.014, 0.35])
# cbar = plt.colorbar(splot, cax=ax2, orientation='vertical', extend='max')
# cbar.set_label('m s$^{-1}$')

# splot.set_visible(False)
plt.tight_layout(pad=3)

figname = "figures/run-average-streamplot_%s_%sm.png" %(expt, zlev)
plt.savefig(figname, dpi=150)
print "......showing"
plt.show()





