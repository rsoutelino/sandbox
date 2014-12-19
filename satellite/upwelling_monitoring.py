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
        self.lon = self.ncfile.variables['lon'][:]
        self.lat = self.ncfile.variables['lat'][:]
        self.lon, self.lat = np.meshgrid(self.lon, self.lat) 
        self.lon = np.ma.masked_where(self.lon > 360, self.lon)
        self.lon = np.ma.masked_where(self.lon < -360, self.lon)
        self.sst = self.ncfile.variables['wtmp'][0,...] - 273.5 # K2C conversion
        
    def plot_image(self, m, meridians, parallels, vmin, vmax, logo_position, lims):
        mlon, mlat = m(self.lon, self.lat)
        fig = plt.figure(facecolor='w', figsize=(15,10))
        ax = fig.add_subplot(111)
        self.ax = ax
        pc = m.pcolormesh(mlon, mlat, self.sst, vmin=vmin, vmax=vmax, 
                     cmap=plt.cm.jet)
        m.fillcontinents()
        # m.bluemarble()
        m.drawcoastlines()
        m.drawstates()
        m.drawmeridians(meridians, labels=[0,0,0,1], dashes=[1,1000])
        m.drawparallels(parallels, labels=[1,0,0,0], dashes=[1,1000])
        ax1 = fig.add_axes([0.17, 0.83, 0.3, 0.02])
        cbar = plt.colorbar(pc, cax=ax1, orientation='horizontal', extend='both')
        self.set_title()
        self.plot_isobaths(ax, lims)
        self.draw_logo(fig, logo_position)

    def plot_isobaths(self, ax, lims):
        topo = sp.loadmat('/home/rsoutelino/misc/etopo/etopo1_atlsul.mat')
        xb, yb, zb = topo['lon'], topo['lat'], topo['topo']
        xb, yb = np.meshgrid(xb, yb)
        xb, yb, zb = romslab.subset(xb, yb, zb, lims)
        cs = ax.contour(xb, yb, -zb, [100, 200, 1000], colors='k', alpha=0.7)
        plt.clabel(cs, fmt='%i')

    def draw_logo(self, fig, logo_position):
        """
        fig : matplotlib figure object
        axes: list with axes position and size
        """              
        ax2 = fig.add_axes(logo_position)
        logofile = Image.open('logo_SR.png')
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
        self.sst = imagelist[0].sst.copy()
        self.avg = np.ones(self.sst.shape)
        for image in imagelist[1:]:
            sst1, sst1_refil = self.sst.copy(), self.sst.copy()
            sst2, sst2_refil = image.sst.copy(), image.sst.copy()
            bad1 = np.where(sst1==False)
            bad2 = np.where(sst2==False)
            sst1_refil[bad1] = sst2[bad1]
            sst2_refil[bad2] = sst1[bad2]

            self.sst = (sst1_refil + sst2_refil) / 2
            # good = np.where(self.sst != False)
            # self.avg[good] = self.avg[good] + 1

        # self.sst = self.sst / self.avg
        self.lon = imagelist[0].lon
        self.lat = imagelist[0].lat
    
    def set_title(self):
        plt.title(u"MSG-SST: Mosaico Diário de %s-%s-%s" %(self.day[:4], 
                     self.day[4:6], self.day[6:8]), fontsize=10)


                     
################################################################################

# lims = (-70, 0, -60, 20)
# lims = (-50, -36, -30, -18)
lims = (-46, -39, -26, -21)
print "loading basemap..............."
m1 = Basemap(projection='cyl', llcrnrlon=lims[0], urcrnrlon=lims[1],
                    llcrnrlat=lims[2], urcrnrlat=lims[3], lat_ts=0, resolution='i')

# METEOSAT SST ############################

pathname = "/home/rsoutelino/ieapm/eumetsat/"
pattern = "201401"
filelist2 = [] # list of total files matching above pattern
imagelist = []
filelist = glob.glob("*0000*.nc")



filelist.sort()

for File in filelist:
    print File
    filelist2.append(File)
    msg = MsgSST(pathname + File, File)
    msg.lon, msg.lat, msg.sst = romslab.subset(msg.lon, msg.lat, msg.sst, lims) 
    imagelist.append(msg)  
    logo_position = [0.134, 0.6, 0.2, 0.2]     
    msg.plot_image(m1, np.arange(-90, 10, 2), np.arange(-70, 30, 2),
                   17, 29, logo_position, lims)
    plt.savefig("upwelling_figs/" + File + ".png", dpi=90)
    plt.close('all')
        
mosaic = Mosaic(imagelist)
mosaic.plot_image(m1, np.arange(-90, 10, 2), np.arange(-70, 30, 2), 18, 29,
                  logo_position, lims)
plt.savefig("upwelling_figs/mosaico_%s.png" %mosaic.day, dpi=90)
plt.close('all')

mdict = {'lon':mosaic.lon, 'lat':mosaic.lat, 'sst':mosaic.sst}
sp.savemat('mosaico_%s.mat' %mosaic.day, mdict)







