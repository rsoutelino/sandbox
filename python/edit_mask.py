#!/usr/bin/env python
######################################################
## Edits ROMS masks using mouse click
## Mar/2012, EDDIES Group - IEAPM
## rsoutelino@gmail.com
######################################################
import wx
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.io as sp
import netCDF4 as nc
from   mpl_toolkits.basemap import Basemap
import romslab

### OPTIONS ##########################################################

# grdname = wx.FileSelector()
grdname = '/home/rsoutelino/myroms/cf_run/exp2/cf2_grd.nc'

### DEFINING FUNCTIONS ###################################

def tellme(s):
    print s
    plt.title(s, fontsize=12)
    plt.show()

### SCRIPT STARTS HERE ########################################################

grd = romslab.RomsGrid(grdname)

### DEFINING MAP PROJECTIONS #################################################

print "\n\n LOADING MAP PROJECTION %s\n\n Should take a while ...\n\n" % ("."*50)

a, b, c, d = grd.corners()
m = Basemap(projection='merc',llcrnrlon=a, urcrnrlon=b,
            llcrnrlat=c, urcrnrlat=d, lat_ts=0, resolution='h')
del a, b, c, d

grd.mlonr, grd.mlatr = m(grd.lonr, grd.latr)

m.pcolormesh(grd.mlonr, grd.mlatr, grd.maskr, vmin=0, vmax=2, cmap=plt.cm.RdBu_r)
m.drawcoastlines()
m.plot(grd.mlonr, grd.mlatr, 'k', alpha=0.5)
m.plot(grd.mlonr.transpose(), grd.mlatr.transpose(), 'k', alpha=0.5)
plt.show()

### Masking phase
c = 'y'

while c == 'y':

    tellme('Pan/Zoom the desired area and go to the prompt:')
    p = raw_input('Continue? [y/n] ')
    tellme('Left click on the lower left corner of the water cell you want to mask:')

    while p != []:
        p = plt.ginput(n=1, timeout=0)
        if p == []:
            break
        line, col = romslab.near2d(grd.mlonr, grd.mlatr, p[0][0], p[0][1])
        grd.maskr[line, col] = 0
        plt.fill( (grd.mlonr[line,col], grd.mlonr[line+1,col], grd.mlonr[line+1,col+1],
                grd.mlonr[line,col+1], grd.mlonr[line,col] ),
                (grd.mlatr[line,col], grd.mlatr[line+1,col], grd.mlatr[line+1,col+1], 
                grd.mlatr[line,col+1], grd.mlatr[line,col] ), 'b' )
        tellme('Keep masking! Just click the scroll buttom when you are done.')

    tellme('Do you wish to mask another area? Go to the prompt. ')
    c = raw_input('Do you wish to mask another area? [y/n] ')

### Unmasking phase
c = 'y'
tellme("Let's start the unmasking phase then. Pan/Zoom in and go to the prompt:")

while c == 'y':
    p = raw_input('Continue? [y/n] ')
    tellme('Time to unmask now. Please click the land cell you wish to unmask:')

    while p != []:
        p = plt.ginput(n=1, timeout=0)
        if p == []:
            break
        line, col = romslab.near2d(grd.mlonr, grd.mlatr, p[0][0], p[0][1])
        grd.maskr[line, col] = 1
        plt.fill( (grd.mlonr[line,col], grd.mlonr[line+1,col], grd.mlonr[line+1,col+1],
                grd.mlonr[line,col+1], grd.mlonr[line,col] ),
                (grd.mlatr[line,col], grd.mlatr[line+1,col], grd.mlatr[line+1,col+1], 
                grd.mlatr[line,col+1], grd.mlatr[line,col] ), 'w' )
        tellme('Keep unmasking! Just click the scroll buttom when you are done.')
        
    tellme('Do you wish to unmask another area? Go to the prompt. ')
    c = raw_input('Do you wish to unmask another area? [y/n] ')



### Creating masku and maskv based on maskr
[grd.masku, grd.maskv, grd.maskp] = romslab.uvp_mask(grd.maskr)

### Writing changes to netcdf file
grd.ncfile.variables['mask_rho'][:] = grd.maskr
grd.ncfile.variables['mask_u'][:] = grd.masku
grd.ncfile.variables['mask_v'][:] = grd.maskv
grd.ncfile.variables['mask_psi'][:] = grd.maskp

grd.ncfile.sync()


















