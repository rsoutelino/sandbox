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
import scipy.io as sp
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near


# PYTHON SCRIPT START ######################################################

# ROMS SETTINGS ==================================================================

expt        = 'phd15'      # case nickname
filetype    = 'avg'        # could be his, avg, rst, ini, clm, ...
days        = np.arange(29,30) # length of the series

nk          = 100           # vertical level [m]
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


print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
pathname = "/home/rsoutelino/myroms/phd_run/"
roms = run_setup(pathname + expt + '_run.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(roms.rundir + roms.run_name + '_grd.nc')

# assigning some variables from grid file
lon   = grdfile.variables['lon_rho'][:]
lat   = grdfile.variables['lat_rho'][:]
h     = grdfile.variables['h'][:]

grdfile  = netCDF4.Dataset(roms.rundir + 'phd13_grd.nc')
h2     = grdfile.variables['h'][:]

fig1 = plt.figure(facecolor='w', figsize=(16,10))
m = Basemap(projection='merc',
    llcrnrlat = lat.min(), urcrnrlat = lat.max(), 
    llcrnrlon = lon.min(), urcrnrlon = lon.max(),
      lat_ts=0, resolution='i')
lonm, latm = m(lon,lat)

plt.subplot(121)
con = m.contourf(lonm, latm, h, 30, cmap=plt.cm.Blues, alpha=0.5)
m.contour(lonm,latm,-h,[-200 -1000],color='k')
#m.plot(np.array([lonm.min(), lonm.min(), lonm.max(), lonm.max(), lonm.min()]),
       #np.array([latm.min(), latm.max(), latm.max(), lonm.min(), lonm.min()]), 'k', linewidth=2)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(-28, -5, 2),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(-44, -30, 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
    
tx, ty = m(-40.4,-11)
plt.text(tx,ty,'Realistic Topography')
#tx, ty = m(-40.4,-11.7)
#plt.text(tx,ty,'Domain')
tx, ty = m(-36.2,-12)
plt.text(tx,ty,'$Lm = 179$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-13)
plt.text(tx,ty,'$Mm = 271$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-14)
plt.text(tx,ty,'$dx = dy = 1/24^\circ$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-15)
plt.text(tx,ty,'$20\ s-levels$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-16)
plt.text(tx,ty,'$hmax = 1500 m$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-17)
plt.text(tx,ty, r'$\frac{\Delta h}{h} = 0.2$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = m(-33,-22)
plt.text(ax, ay, "a", ha="center", va="center", size=12, bbox=bbox_props)

ax = fig1.add_axes([0.15, 0.4, 0.01, 0.4])
cbar = plt.colorbar(con,cax=ax,orientation='vertical',alpha=0.2)
cbar.set_label('Depth [m]')


plt.subplot(122)
con = m.contourf(lonm, latm, h2, 30, cmap=plt.cm.Blues, alpha=0.5)
m.contour(lonm,latm,-h,[-200 -1000],color='k')
#m.plot(np.array([lonm.min(), lonm.min(), lonm.max(), lonm.max(), lonm.min()]),
       #np.array([latm.min(), latm.max(), latm.max(), lonm.min(), lonm.min()]), 'k', linewidth=2)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(-28, -5, 2),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(-44, -30, 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
    
tx, ty = m(-40.4,-11)
plt.text(tx,ty,'Flat Bottom')
#tx, ty = m(-40.4,-11.7)
#plt.text(tx,ty,'Domain')
tx, ty = m(-36.2,-12)
plt.text(tx,ty,'$Lm = 179$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-13)
plt.text(tx,ty,'$Mm = 271$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-14)
plt.text(tx,ty,'$dx = dy = 1/24^\circ$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-15)
plt.text(tx,ty,'$20\ s-levels$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-16)
plt.text(tx,ty,'$hmax = 1500 m$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})
tx, ty = m(-36.2,-17)
plt.text(tx,ty, r'$\frac{\Delta h}{h} = 0.2$', {'fontsize' : 15, 'color' : 'w', 'fontweight' : 'bold'})

bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
ax,ay = m(-33,-22)
plt.text(ax, ay, "b", ha="center", va="center", size=12, bbox=bbox_props)

ax = fig1.add_axes([0.57, 0.4, 0.01, 0.4])
cbar = plt.colorbar(con,cax=ax,orientation='vertical',alpha=0.2)
cbar.set_label('Depth [m]')

plt.show()
plt.savefig('figures/model_grid.pdf') 

stop


del m

m = Basemap(projection='merc',
    llcrnrlat = -35, urcrnrlat = 0, 
    llcrnrlon = -54, urcrnrlon = 20,
      lat_ts=0, resolution='l')
lonm, latm = m(lon,lat)

etopo = sp.loadmat('/home/rsoutelino/rsoutelino/phd/data/etopo2_leste.mat')
xb = etopo.pop('xb')
yb = etopo.pop('yb')
zb = etopo.pop('zb')
mxb,myb = m(xb,yb)

#plt.figure(facecolor='w')
#m.pcolormesh(mxb,myb,zb,cmap=plt.cm.PuBu)
##m.contour(mxb,myb,zb,[-200 -1000],color='k')
##m.contourf(lonm, latm, h, 30, cmap=plt.cm.Blues, alpha=0.2)
#plt.plot(np.array([lonm.min(), lonm.min(), lonm.max(), lonm.max(), lonm.min()]),
       #np.array([latm.min(), latm.max(), latm.max(), lonm.min(), lonm.min()]), 'k', linewidth=2)
#m.drawcoastlines(zorder=5)
#m.drawcountries(zorder=4)
#m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
#m.drawparallels(np.arange(-40, 0, 5),
    #labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6, color='0.6')
#m.drawmeridians(np.arange(-60, 20, 10),
    #labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7, color='0.6')
#plt.show()
#plt.savefig('figures/model_grid_large.pdf') 
    
stop
m = Basemap(projection='merc',
    llcrnrlat = -23, urcrnrlat = -12, 
    llcrnrlon = -42.5, urcrnrlon = -30,
      lat_ts=0, resolution='i')

mxb,myb = m(xb,yb)


fig3 = plt.figure(facecolor='w')
con = m.pcolormesh(mxb, myb, -zb, vmin=0, vmax=5000, cmap=plt.cm.Spectral)
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
m.drawparallels(np.arange(-28, -5, 2),
    labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
m.drawmeridians(np.arange(-44, -30, 2),
    labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)
    
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.5)
a = plt.ginput(n=3)
plt.text(a[0][0], a[0][1], "Royal-Charlotte Bank (RCB)", ha="left", va="center", size=12, bbox=bbox_props)
plt.text(a[1][0], a[1][1], "Abrolhos Bank (AB)", ha="left", va="center", size=12, bbox=bbox_props)
plt.text(a[2][0], a[2][1], "Vitoria-Trindade Ridge (VTR)", ha="left", va="center", size=12, bbox=bbox_props)

ax = fig3.add_axes([0.24, 0.35, 0.015, 0.5])
cbar = plt.colorbar(con,cax=ax,orientation='vertical',alpha=0.2)
cbar.set_label('Depth [m]')

plt.savefig('figures/topog.pdf')

plt.show()



