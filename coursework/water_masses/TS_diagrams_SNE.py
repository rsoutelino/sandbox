#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# TS diagrams in S, N and E boundaries of BC origin domain
#
# Rafael Soutelino - rsoutelino@gmail.com
#
#
# Last modification: Apr, 2011
#####################################################################

print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 

# IMPORTING MODULES #################################################
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# import datetime as dt
import netCDF4 as nc
import scipy.io as sp
import seawater.csiro as sw

#####################################################################

# base to plot isopicnals
Sg = np.arange(33, 38, 0.25)
Tg = np.arange(0, 30)
Sg, Tg = np.meshgrid(Sg, Tg)
Dg = sw.dens0(Sg, Tg) - 1000

at = 20
acas = (20, 8.72)
aia = (8.72, 3.31)
apan = (3.31, 2.04)

plt.figure(1,figsize=(16,10),facecolor='w')

# loading data
dataset = sp.loadmat('oe2_SNE.mat')
tempSoe2 = dataset['tempSoe2']; saltSoe2 = dataset['saltSoe2'];
tempNoe2 = dataset['tempNoe2']; saltNoe2 = dataset['saltNoe2'];
tempEoe2 = dataset['tempEoe2']; saltEoe2 = dataset['saltEoe2'];

saltSoe2 = saltSoe2.ravel(); tempSoe2 = tempSoe2.ravel();
f1 = np.where(tempSoe2 > at)
f2 = np.where( (tempSoe2 < acas[0]) & (tempSoe2 > acas[1]) )
f3 = np.where( (tempSoe2 < aia[0])  & (tempSoe2 > aia[1])  )
f4 = np.where( (tempSoe2 < apan[0]) & (tempSoe2 > apan[1]) )

p1 = plt.subplot(331)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltSoe2[f1], tempSoe2[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltSoe2[f2], tempSoe2[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltSoe2[f3], tempSoe2[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltSoe2[f4], tempSoe2[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])

saltNoe2 = saltNoe2.ravel(); tempNoe2 = tempNoe2.ravel();
f1 = np.where(tempNoe2 > at);
f2 = np.where( (tempNoe2 < acas[0]) & (tempNoe2 > acas[1]) )
f3 = np.where( (tempNoe2 < aia[0])  & (tempNoe2 > aia[1])  )
f4 = np.where( (tempNoe2 < apan[0]) & (tempNoe2 > apan[1]) )

p2 = plt.subplot(332)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltNoe2[f1], tempNoe2[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltNoe2[f2], tempNoe2[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltNoe2[f3], tempNoe2[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltNoe2[f4], tempNoe2[f4], 'm.', markersize=1, alpha=0.1)
p2.set_title('OE II')
plt.axis([33.5, 37.5, 0, 30])

saltEoe2 = saltEoe2.ravel(); tempEoe2 = tempEoe2.ravel();
f1 = np.where(tempEoe2 > at);
f2 = np.where( (tempEoe2 < acas[0]) & (tempEoe2 > acas[1]) )
f3 = np.where( (tempEoe2 < aia[0])  & (tempEoe2 > aia[1])  )
f4 = np.where( (tempEoe2 < apan[0]) & (tempEoe2 > apan[1]) )

p3 = plt.subplot(333)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltEoe2[f1], tempEoe2[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltEoe2[f2], tempEoe2[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltEoe2[f3], tempEoe2[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltEoe2[f4], tempEoe2[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])


# loading data
dataset = sp.loadmat('woa2009_mar_SNE.mat')
tempSwoa = dataset['tempSwoa']; saltSwoa = dataset['saltSwoa'];
tempNwoa = dataset['tempNwoa']; saltNwoa = dataset['saltNwoa'];
tempEwoa = dataset['tempEwoa']; saltEwoa = dataset['saltEwoa'];

saltSwoa = saltSwoa.ravel(); tempSwoa = tempSwoa.ravel();
f1 = np.where(tempSwoa > at);
f2 = np.where( (tempSwoa < acas[0]) & (tempSwoa > acas[1]) )
f3 = np.where( (tempSwoa < aia[0])  & (tempSwoa > aia[1])  )
f4 = np.where( (tempSwoa < apan[0]) & (tempSwoa > apan[1]) )

p4 = plt.subplot(334)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltSwoa[f1], tempSwoa[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltSwoa[f2], tempSwoa[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltSwoa[f3], tempSwoa[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltSwoa[f4], tempSwoa[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])

saltNwoa = saltNwoa.ravel(); tempNwoa = tempNwoa.ravel();
f1 = np.where(tempNwoa > at);
f2 = np.where( (tempNwoa < acas[0]) & (tempNwoa > acas[1]) )
f3 = np.where( (tempNwoa < aia[0])  & (tempNwoa > aia[1])  )
f4 = np.where( (tempNwoa < apan[0]) & (tempNwoa > apan[1]) )

p5 = plt.subplot(335)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltNwoa[f1], tempNwoa[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltNwoa[f2], tempNwoa[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltNwoa[f3], tempNwoa[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltNwoa[f4], tempNwoa[f4], 'm.', markersize=1, alpha=0.1)
p5.set_title('WOA2009 - March')
plt.axis([33.5, 37.5, 0, 30])

saltEwoa = saltEwoa.ravel(); tempEwoa = tempEwoa.ravel();
f1 = np.where(tempEwoa > at);
f2 = np.where( (tempEwoa < acas[0]) & (tempEwoa > acas[1]) )
f3 = np.where( (tempEwoa < aia[0])  & (tempEwoa > aia[1])  )
f4 = np.where( (tempEwoa < apan[0]) & (tempEwoa > apan[1]) )

p6 = plt.subplot(336)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltEwoa[f1], tempEwoa[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltEwoa[f2], tempEwoa[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltEwoa[f3], tempEwoa[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltEwoa[f4], tempEwoa[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])

# loading data
dataset = sp.loadmat('occam2003_mar_SNE.mat')
tempSoccam = dataset['tempSoccam']; saltSoccam = dataset['saltSoccam'];
tempNoccam = dataset['tempNoccam']; saltNoccam = dataset['saltNoccam'];
tempEoccam = dataset['tempEoccam']; saltEoccam = dataset['saltEoccam'];

saltSoccam = saltSoccam.ravel(); tempSoccam = tempSoccam.ravel();
f1 = np.where(tempSoccam > at);
f2 = np.where( (tempSoccam < acas[0]) & (tempSoccam > acas[1]) )
f3 = np.where( (tempSoccam < aia[0])  & (tempSoccam > aia[1])  )
f4 = np.where( (tempSoccam < apan[0]) & (tempSoccam > apan[1]) )

p7 = plt.subplot(337)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltSoccam[f1], tempSoccam[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltSoccam[f2], tempSoccam[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltSoccam[f3], tempSoccam[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltSoccam[f4], tempSoccam[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])
p7.set_xlabel('South')

saltNoccam = saltNoccam.ravel(); tempNoccam = tempNoccam.ravel();
f1 = np.where(tempNoccam > at);
f2 = np.where( (tempNoccam < acas[0]) & (tempNoccam > acas[1]) )
f3 = np.where( (tempNoccam < aia[0])  & (tempNoccam > aia[1])  )
f4 = np.where( (tempNoccam < apan[0]) & (tempNoccam > apan[1]) )

p8 = plt.subplot(338)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltNoccam[f1], tempNoccam[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltNoccam[f2], tempNoccam[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltNoccam[f3], tempNoccam[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltNoccam[f4], tempNoccam[f4], 'm.', markersize=1, alpha=0.1)
p8.set_title('OCCAM2003 - March')
plt.axis([33.5, 37.5, 0, 30])
p8.set_xlabel('North')

saltEoccam = saltEoccam.ravel(); tempEoccam = tempEoccam.ravel();
f1 = np.where(tempEoccam > at);
f2 = np.where( (tempEoccam < acas[0]) & (tempEoccam > acas[1]) )
f3 = np.where( (tempEoccam < aia[0])  & (tempEoccam > aia[1])  )
f4 = np.where( (tempEoccam < apan[0]) & (tempEoccam > apan[1]) )

p9 = plt.subplot(339)
c=plt.contour(Sg, Tg, Dg, range(20, 40), cmap=None, colors='0.8')
c.clabel(fmt='%1.0f',fontsize=10)
plt.plot(saltEoccam[f1], tempEoccam[f1], 'r.', markersize=1, alpha=0.1)
plt.plot(saltEoccam[f2], tempEoccam[f2], 'g.', markersize=1, alpha=0.1)
plt.plot(saltEoccam[f3], tempEoccam[f3], 'b.', markersize=1, alpha=0.1)
plt.plot(saltEoccam[f4], tempEoccam[f4], 'm.', markersize=1, alpha=0.1)
plt.axis([33.5, 37.5, 0, 30])
p9.set_xlabel('East')

plt.show()


