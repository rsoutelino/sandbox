#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# map.py
#
# purpose:  quick-and-dirty map using matplotlib's basemap toolkit
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  02-Sep-2010
# modified: Tue 02 Nov 2010 04:33:57 PM EDT
#
# obs:  map lonlat.dat -87.5 -22.5 -59.5 14.5
#

import sys
import os
import matplotlib
import numpy

# matplotlib.use('Agg')

""" import basemap and pylab """
try:
  from mpl_toolkits.basemap import Basemap
except ImportError:
    print "Basemap not found."

try:
  from pylab import *
except ImportError:
    print "Matplotlib not found."

__version__ = "0.3.2"
# parse args (ugly!)
if len(sys.argv) == 7:
    lonmin  = float(sys.argv[3])
    lonmax  = float(sys.argv[4])
    latmin  = float(sys.argv[5])
    latmax  = float(sys.argv[6])
    lonlat1 = np.loadtxt(sys.argv[1])
    lonlat2 = np.loadtxt(sys.argv[2])
    lonpt1  = lonlat1[:,0];
    lonpt2  = lonlat2[:,0];
    latpt1  = lonlat1[:,1];
    latpt2  = lonlat2[:,1];
elif len(sys.argv) == 6:
    lonmin  = float(sys.argv[2])
    lonmax  = float(sys.argv[3])
    latmin  = float(sys.argv[4])
    latmax  = float(sys.argv[5])
    lonlat1 = np.loadtxt(sys.argv[1])
    lonpt1  = lonlat1[:,0]
    latpt1  = lonlat1[:,1]
elif len(sys.argv) == 5:
    lonmin = float(sys.argv[1])
    lonmax = float(sys.argv[2])
    latmin = float(sys.argv[3])
    latmax = float(sys.argv[4])
else:
    sys.exit("""\nMust provide at least 4 arguments!!!
             A 2 column file with longitue and latitute of the points (optional)
             and the lon lat limits for the map.
             e.g.: map lonlat.dat -87.5 -22.5 -59.5 14.5\n""")

m = Basemap( projection = 'merc',
             llcrnrlat = latmin, urcrnrlat = latmax,
             llcrnrlon = lonmin, urcrnrlon = lonmax,
             resolution = 'l' ) # low coastlines resolution

""" create figure and associate axes with Basemap """
fig  = figure()
ax   = fig.add_subplot(111)
m.ax = ax

#shp_info = m.readshapefile('/usr/share/basemap/brazil','states',drawbounds=True)
#m.drawstates()
m.drawcountries(linewidth=1.2)
m.drawcoastlines(linewidth=0.6)
m.bluemarble()
m.drawparallels(np.arange(-25,-12,1),labels=[1, 1, 0, 0],
    dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-44,-32,1),labels=[1, 1, 1, 1],dashes=[1,1000],
    zorder=7)

""" add a non bluearth alternative """
#m.drawcoastlines()
#m.drawlsmask(land_color='grey',ocean_color='None', lakes=True)

""" plot points """
msize=5 #size of the points
if len(sys.argv) == 6:
    xc1,yc1 = m(lonpt1,latpt1)
    #m.plot(xc1,yc1,'ro', alpha = 0.5, markersize=2)
    m.plot(xc1, yc1, color='red', marker='o', markersize=msize, linewidth=0)
    if lonpt1.min() < lonmin:
        print "\n some points are out of bounds (longitude min is = %f and limit was %f) \n " % (lonpt1.min(),lonmin)
    if lonpt1.max() > lonmax:
        print "\n some points are out of bounds (longitude max is = %f and limit was %f) \n " % (lonpt1.max(),lonmax)
    if latpt1.min() < latmin:
        print "\n some points are out of bounds (latitude min is = %f and limit was %f) \n "  % (latpt1.min(),latmin)
    if latpt1.max() > latmax:
        print "\n some points are out of bounds (latitude max is = %f and limit was %f) \n "  % (latpt1.max(),latmax)
    print "WARNING: \n only 1 lon/lat file to plot \n"
elif len(sys.argv) == 7:
    xc1,yc1 = m(lonpt1,latpt1)
    xc2,yc2 = m(lonpt2,latpt2)
    m.plot(xc1, yc1, color='red', marker='o', markersize=msize, linewidth=0) 
    m.plot(xc2, yc2, color='yellow', marker='o', markersize=msize, linewidth=0)
    if lonpt1.min() < lonmin:
        print "\n some points are out of bounds in the first file (longitude min is = %f and limit was %f) \n " % (lonpt1.min(), lonmin)
    if lonpt1.max() > lonmax:
        print "\n some points are out of bounds in the first file (longitude max is = %f and limit was %f) \n " % (lonpt1.max(), lonmax)
    if latpt1.min() < latmin:
        print "\n some points are out of bounds in the first file (latitude min is = %f and limit was %f) \n "  % (latpt1.min(), latmin)
    if latpt1.max() > latmax:
        print "\n some points are out of bounds in the first file (latitude max is = %f and limit was %f) \n "  % (latpt1.max(), latmax)
    if lonpt2.min() < lonmin:
        print "\n some points are out of bounds in the first file (longitude min is = %f and limit was %f) \n " % (lonpt2.min(), lonmin)
    if lonpt2.max() > lonmax:
        print "\n some points are out of bounds in the first file (longitude max is = %f and limit was %f) \n " % (lonpt2.max(), lonmax)
    if latpt2.min() < latmin:
        print "\n some points are out of bounds in the first file (latitude min is = %f and limit was %f) \n "  % (latpt2.min(), latmin)
    if latpt2.max() > latmax:
        print "\n some points are out of bounds in the first file (latitude max is = %f and limit was %f) \n "  % (latpt2.max(), latmax)
    print "WARNING: \n two lon/lat files to plot \n"
else:
    print "WARNING: no lon/lat files to plot \n"

#show()

""" trim image """
import StringIO, Image
imgdata = StringIO.StringIO()
fig.savefig(imgdata, dpi=300, format='png')
imgdata.seek(0)  # rewind the data
im = Image.open(imgdata)

def trim(im, border):
  from PIL import ImageChops
  bg = Image.new(im.mode, im.size, border)
  diff = ImageChops.difference(im, bg)
  bbox = diff.getbbox()
  if bbox:
      return im.crop(bbox)
  else:
      # found no content
      raise ValueError("cannot trim; image was empty")

im = trim(im,'white')
im.show()
