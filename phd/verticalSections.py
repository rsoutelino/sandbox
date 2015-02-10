#!/usr/bin/env python
# -*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.io as sp
from mpl_toolkits.basemap import Basemap
import datetime as dt
import numpy as np
import matplotlib
import romslab

matplotlib.rcParams.update({'font.size': 10})
####################################################################

# def plotSection(Roms):
#     fig = plt.

####################################################################

pathname = "/home/rsoutelino/myroms/phd_run/"
outputList = ["phd15_avg.nc", "phd18_avg.nc", "phd16_avg.nc"]

# days = [1, 10, 20, 30]
days = [10, 30, 60, 90]
secStart = [ (-39, -12), (-39, -15), (-39, -17), (-39, -19) ]
secStart = [ (-40, -12), (-40, -15), (-40, -17), (-40, -19) ]
secEnd   = [ (-36, -12), (-36, -15), (-36, -17), (-36, -19) ]


s = 0.05
contours = np.arange(-0.3, 0.3 + s, s)
zlim = (-1500, 0)
cm = plt.cm.RdBu
dayLocator = (secStart[0][0] + 0.2, zlim[0] + 200) 
latLocator = (secStart[0][0] + 0.2, zlim[0] + 100)

i = 2
print "\nInstanciating PlotROMS class ........\n"
Roms = romslab.PlotROMS( pathname + outputList[i] )

fig, axarr = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(10,12))
fig.set_facecolor('w')

for l in range(len(days)):
    for y in range(len(secStart)):
        Roms.vslice(secStart[y], secEnd[y], contours, 
                zlim, 'v', cmap=plt.cm.RdBu, l=days[l]-1)
        Roms.vVslice = np.ma.masked_where(np.abs(Roms.vVslice) > 1, Roms.vVslice)
        ax = axarr[y, l]
        ax.set_axis_bgcolor('0.7')
        con = ax.contourf(Roms.xVslice, Roms.zVslice, Roms.vVslice, contours,
                    cmap=cm, extend='both')
        for label in ax.xaxis.get_ticklabels():
            label.set_rotation(45)
        ax.text(dayLocator[0], dayLocator[1], "Day %s" %days[l] )
        ax.text(latLocator[0], latLocator[1], r"%s S" %np.abs(secStart[y][1]) )
        if y==3: ax.set_xlabel('Longitude')
        if l==0: ax.set_ylabel('Z [m]')

axcbar = fig.add_axes([0.1, 0.82, 0.015, 0.15])
cbar = plt.colorbar(con, cax=axcbar)
cbar.set_label('m s$^{-1}$')

plt.tight_layout(pad=0.5)
plt.show()
plt.savefig('figures/sec_vertical_%s.pdf' %outputList[i][:5])




