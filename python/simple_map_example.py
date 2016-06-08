# -*- coding: utf-8 -*-
######################################################
## define title!!!!!!!!!!!  date, authorship        ##
######################################################
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

######################################################

xmin, xmax, ymin, ymax = (-41, -35, -25, -20)

m = Basemap(projection='cyl', llcrnrlat=ymin, urcrnrlat=ymax,
    llcrnrlon=xmin, urcrnrlon=xmax, lat_ts=0, resolution='i')

data = nc.Dataset('/path/nomedoarquivo.nc')
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
prop = data.variables['prop'][:]

#
# processamento .......
#
#

print lon.shape, lat.shape, prop.shape

# plotting map:
plt.figure()
m.pcolormesh(lon, lat, prop, vmin=xxx, vmax=xxx, cmap=plt.cm.RdBu)
plt.colorbar()
m.fillcontinents()
m.drawcoastlines()
m.drawstates()
m.drawmeridians(meridians, labels=[0,0,0,1], dashes=[1,1000])
m.drawparallels(parallels, labels=[1,0,0,0], dashes=[1,1000])
plt.savefig("filename.png")