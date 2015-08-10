#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  plot_roms.py -- Visualization of ROMS fields
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Jul, 2013
###################################################################################

import numpy as np
import sys
import scipy.io as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import seawater.csiro as sw
import netCDF4 as nc
from matplotlib.patches import Polygon

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near
from cookb_signalsmooth import smooth




def compute_baroclinic_conversion(expt, eddy, zlev):
    avg  = sp.loadmat("eddies_subsets/%s_%s_subset_%sm.mat" %(expt, eddy, zlev))
    lon  = avg['lon']
    lat  = avg['lat']
    u    = avg['u']
    v    = avg['v']
    t    = avg['temp']

    u = np.ma.masked_where(np.abs(u) > 2, u)
    v = np.ma.masked_where(np.abs(v) > 2, v)
    t = np.ma.masked_where(np.abs(v) > 30, t)


    lm, jm, im = u.shape

    # computing averages and perturbations ####################################
    ubar, vbar, tbar = u.mean(axis=0), v.mean(axis=0), t.mean(axis=0)
    unot, vnot, tnot = u*0, v*0, t*0

    for l in range(lm):
        unot[l,...] = u[l,...] - ubar
        vnot[l,...] = v[l,...] - vbar
        tnot[l,...] = t[l,...] - tbar


    # baroclinic energy  conversions ##########################################
    bec = []
    g = 9.8

    theta_z = tbar.mean()
    alpha = sw.alpha(35, theta_z, zlev)
    gradT = np.gradient(tbar)
    dTdx  = gradT[1] / (0.05 * 111000)
    dTdy  = gradT[0] / (0.05 * 111000)

    count = 0

    for l in range(lm):
        count += 1
        bc =  ( (g*alpha)/theta_z ) * ( unot[l]*tnot[l]*dTdx + vnot[l]*tnot[l]*dTdy )
        bec.append( bc.mean() )

    bec = np.array(bec)

    bec = smooth(bec, window_len=10, window='hanning')

    return bec




def make_subplot(sbplt, zlev, bec, time, color, title='no', xlabel='no', eddy='IE'):
    p = plt.subplot(sbplt)

    if title == 'yes':
        title = "Horizontally-averaged 5-years Baroclinic Conversion Rate\
 for %s Eddy" %(eddy)
        plt.title(title, fontsize=10, fontweight='bold')

    plt.plot(time, bec, color, label=eddy, linewidth=2)
    f = np.where(bec < 0)
    plt.plot(time[f], bec[f], 'w--', linewidth=3)
    becFill  = bec.copy()
    becFill[ np.where(bec < 0) ] = 0
    verts = [(time.min(), 0)] + zip(time, becFill) + [(time.max(), 0)]
    poly = Polygon(verts, facecolor=color, edgecolor=color, alpha=0.3)
    p.add_patch(poly)



    plt.axis([5, 10, -5e-12, 5e-12])
    plt.legend()
    plt.ylabel('m$^2$ s$^{-3}$', fontweight='bold')
    
    if xlabel == 'yes':
        plt.xlabel('Years', fontweight='bold')

    plt.grid()
    plt.show()



def plot_days_annotations():
    # plotting selected snapshots time
    days = [1886, 2555, 2920, 3406]
    times = [5.16, 7.0, 8.0, 9.33]
    strlist = ["Year 5\nDay 60", "Year 7\nDay 1", "Year 8\nDay 1", "Year 9\nDay 120"]
    for x, date in zip(times, strlist):
        head_length = 7e-13
        width = 0.1
        y = -2e-12
        dx = 0.0
        dy = 2e-12 - head_length
        # if x == 7.0:
        #     plt.arrow(x, -y, dx, -dy, color='k', alpha=0.4, shape='full', width=width, 
        #                               head_width=0.2, head_length=7e-13, zorder=11)
        #     plt.text(x-width, y*-1 + 4e-13 , date, fontsize=9, fontweight='bold', zorder=12)
        # else:
        plt.arrow(x, y, dx, dy, color='0.7', shape='full', width=width, 
                                head_width=0.2, head_length=7e-13, zorder=11)
        plt.text(x-width, y-1e-12, date, fontsize=9, fontweight='bold', zorder=12)



# PYTHON SCRIPT START ######################################################

plt.close('all')

expt = 'phd16'

zlevs = [50, 100, 400]

plt.figure(facecolor='w', figsize=(10,14))

eddy = 'IE'
bec = compute_baroclinic_conversion(expt, eddy, zlevs[0])
time = np.linspace(5, 10, bec.size)
make_subplot(sbplt=311, zlev=zlevs[0], bec=bec, time=time, color='b', title='no', eddy=eddy)
if expt == 'phd15':
    plot_days_annotations()


eddy = 'RCE'
bec = compute_baroclinic_conversion(expt, eddy, zlevs[0])
make_subplot(sbplt=312, zlev=zlevs[1], bec=bec, time=time, color='r', eddy=eddy)
if expt == 'phd15':
    plot_days_annotations()

eddy = 'AbE'
bec = compute_baroclinic_conversion(expt, eddy, zlevs[0])
make_subplot(sbplt=313, zlev=zlevs[2], bec=bec, time=time, color='g', xlabel='yes', eddy='AE')
if expt == 'phd15':
    plot_days_annotations()


figname = "energy_conversion/bc_energy_conversion_%s.pdf" %(expt)
plt.savefig(figname)








