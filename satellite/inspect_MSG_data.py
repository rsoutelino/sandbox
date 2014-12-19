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
import glob
import Image

################################################################################

class MsgSST(object):
    def __init__(self, fullpath, filename):
        self.filename = filename
        self.ncfile = nc.Dataset(fullpath)
        self.time = self.ncfile.variables['time'][:]
        self.lon = self.ncfile.variables['lon'][:] - 360
        self.lat = self.ncfile.variables['lat'][:]
        self.lon, self.lat = np.meshgrid(self.lon, self.lat) 
        self.sst = self.ncfile.variables['wtmp'][0,...] - 273.5 # K2C conversion
        
    def plot_image(self, m, meridians, parallels, vmin, vmax, logo_position):
        mlon, mlat = m(self.lon, self.lat)
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        m.pcolormesh(mlon, mlat, self.sst, vmin=vmin, vmax=vmax)
        # m.fillcontinents()
        m.bluemarble()
        m.drawcoastlines()
        m.drawmeridians(meridians, labels=[0,0,0,1], dashes=[1,1000])
        m.drawparallels(parallels, labels=[1,0,0,0], dashes=[1,1000])
        plt.colorbar()
        self.set_title()
        self.draw_logo(fig, logo_position)
                  
    def draw_logo(self, fig, logo_position):
        """
        fig : matplotlib figure object
        axes: list with axes position and size
        """              
        ax2 = fig.add_axes(logo_position)
        logofile = Image.open('logo-left.png')
        logo = np.asarray(logofile)
        im = plt.imshow(logo)
        ax2.set_axis_off()
    
    def set_title(self):
        plt.title("MSG-SST: %s-%s-%s, %s:%s Z" %(self.filename[30:34], 
                  self.filename[34:36], self.filename[36:38], 
                  self.filename[38:40], self.filename[40:42] ), fontsize=10)
        
        
class Mosaic(MsgSST):
    def __init__(self, imagelist):
        self.imagelist = imagelist
        self.day = str(int(self.imagelist[0].time[:]))
        self.compute_mosaic()
        
    def compute_mosaic(self):
        sst = imagelist[0].sst
        for image in imagelist:
            f = np.where(sst==False)
            sst[f[0], f[1]] = image.sst[f[0], f[1]]  
        self.sst = sst
        self.lon = imagelist[0].lon
        self.lat = imagelist[0].lat
    
    def set_title(self):
        plt.title(u"MSG-SST: Mosaico Diário de %s-%s-%s" %(self.day[:4], 
                     self.day[4:6], self.day[6:8]), fontsize=10)
                     
################################################################################

# lims = (-70, 0, -60, 20)
lims2 = (-46, -39, -26, -20)
lims = lims2
m1 = Basemap(projection='cyl', llcrnrlon=lims[0], urcrnrlon=lims[1],
                    llcrnrlat=lims[2], urcrnrlat=lims[3], lat_ts=0, resolution='l')
m2 = Basemap(projection='cyl', llcrnrlon=lims2[0], urcrnrlon=lims2[1],
                    llcrnrlat=lims2[2], urcrnrlat=lims2[3], lat_ts=0, resolution='i')

# METEOSAT SST ############################

pathname = "/home/rsoutelino/ieapm/eumetsat/"
pattern = "20120918"
filelist2 = [] # list of total files matching above pattern
imagelist = []
filelist = glob.glob("*.nc")
filelist.sort()

for File in filelist:
    if pattern in File and "grb" not in File:
        print File
        filelist2.append(File)
        msg = MsgSST(pathname + File, File)
        msg.lon, msg.lat, msg.sst = romslab.subset(msg.lon, msg.lat, msg.sst, lims2) 
        imagelist.append(msg)       
        # Broad Area ++++++++++++++++++++++++
        msg.plot_image(m1, np.arange(-90, 10, 10), np.arange(-70, 30, 10),
                       -2, 30, [0.6, 0.12, 0.15, 0.15])
        plt.savefig("figures/" + File + ".png")
        plt.close('all')
        # Zoom into Cabo Frio +++++++++++++++
        msg.plot_image(m2, np.arange(-90, 10, 2), np.arange(-70, 30, 2), 16, 26,
                       [0.13, 0.62, 0.15, 0.15])
        plt.savefig("figures/zoom_" + File + ".png")
        plt.close('all')
        
        
mosaic = Mosaic(imagelist)
mosaic.plot_image(m1, np.arange(-90, 10, 10), np.arange(-70, 30, 10), -2, 30,
                    [0.6, 0.12, 0.15, 0.15])
plt.savefig("figures/mosaico_%s.png" %mosaic.day)
mosaic.plot_image(m2, np.arange(-90, 10, 2), np.arange(-70, 30, 2), 16, 26,
                  [0.13, 0.62, 0.15, 0.15])
plt.savefig("figures/mosaico_zoom_%s.png" %mosaic.day)
plt.close('all')

mdict = {'lon':mosaic.lon, 'lat':mosaic.lat, 'sst':mosaic.sst}
sp.savemat('mosaico_%s.mat' %mosaic.day, mdict)







