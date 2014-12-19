#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import seawater.csiro as sw
import scipy.io as sp
import netCDF4 as nc
import romslab

################################################################################

def flatten(nested):
    for sublist in nested:
        for element in sublist:
            yield element

pattern = "S-OSI_-FRA_-MTOP-MGRSST_FIELD"

filelist = os.listdir("data")

lims = (-70, 0, -60, 20)
m = Basemap(projection='cyl', llcrnrlon=lims[0], urcrnrlon=lims[1],
                    llcrnrlat=lims[2], urcrnrlat=lims[3], lat_ts=0, resolution='l')

LON, LAT, SST = [], [], []
filelist2 = [] # list of total files matching above pattern
filelist3 = [] # list of files inside desired area  

for File in filelist:
    if pattern in File:
        data = nc.Dataset("data/" + File)
        lon, lat = data.variables['lon'][::5], data.variables['lat'][::5]
        sst = data.variables['sea_surface_temperature'][::5]
        lon, lat, sst = romslab.subset(lon, lat, sst, lims)       
        if lon != None:
            LON.append( list(lon.ravel()) )
            LAT.append( list(lat.ravel()) )
            SST.append( list(sst.ravel()) )
            print File + "  --> OK!"
            filelist2.append(File)
            filelist3.append(File)         
        else:
            print File + "  --> Out of Range!"
            filelist2.append(File)
     
LON = np.array( list(flatten(LON)) )
LAT = np.array( list(flatten(LAT)) )
SST = np.array( list(flatten(SST)) ) - 273.5

mlon, mlat = m(LON, LAT)

plt.figure(facecolor='w')
m.scatter(mlon[::10], mlat[::10], c=SST[::10], s=10, edgecolors='w',
             linewidth=0.0000001, vmin=-2, vmax=30)
m.fillcontinents()
m.drawcoastlines()
plt.colorbar()
plt.show()









