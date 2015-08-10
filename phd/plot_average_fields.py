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
import netCDF4

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near

plt.close('all')
# PYTHON SCRIPT START ######################################################

# Input parameters
expt = 'phd16'
zlev = 2000
lims = [-41.0, -35, -21, -13]
d    = 3  # subsampling
s    = 5  # quiver size
sc   = 0.5 # quiver scale
spacer  = 0.5 # spacer between scale arrow and scale text 
 
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


plt.figure(figsize=(10,10), facecolor='w')
m.quiver(lon[::d,::d], lat[::d,::d], u[::d,::d]*s, v[::d,::d]*s, scale=10)
m.fillcontinents()
m.drawcoastlines()
m.drawmeridians(np.arange(lims[0], lims[1], 2),
     labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
m.drawparallels(np.arange(lims[2], lims[3], 1),
     labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
# scale
ax, ay = m(lims[0]+spacer, lims[3]-spacer)
m.quiver(ax, ay, sc, 0, scale=10, zorder=9)
plt.text(ax, ay-spacer, "%s m/s" %(sc), color='k',fontsize=10,fontweight='bold')
plt.title("Run-averaged velocities: %s m" %zlev)
figname = "figures/run-average-vels_%s_%sm.png" %(expt, zlev)
plt.savefig(figname)
plt.show()
