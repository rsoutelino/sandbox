#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#	This script has some configurations parameters to set before
#	plotting ROMS fields. 
#
#	Rafael Soutelino - rsoutelino@gmail.com
# 	Last Modification: Aug, 2010
#   Still only working for ROMS OUTPUTS, not INPUTS
#################################################################
import matplotlib.pyplot as plt
import sys

# let's start with the most changed ones:

# ROMS SETTINGS ==================================================================

expt        = 'phd15'      # case nickname
filetype    = 'avg-1'       # could be his, avg, rst, ini, clm, ...
PLOT        = 1           # vel-temp = 1, zonal slice = 2, meridional slice = 3
			 # vbar-zeta = 4, zonal temp = 5, meridional temp = 6
l           = 12         # time in model output [time-step]
movie       = 1          # make movie or snapshot


days = range(360, 500, 30) 



# PLOT = 1 =======================================================================

if PLOT == 1:
    nk          = 50           # vertical level [m]
    lonplot     = (-40, -36)  # longitude limits for plotting
    latplot     = (-19.5, -14)  # latitude limits for plotting
    # lonplot     = (-41, -34)  # longitude limits for plotting
    # latplot     = (-22.5, -12)  # latitude limits for plotting
    #lonplot     = (-41, -32)  # longitude limits for plotting
    #latplot     = (-23, -10)  # latitude limits for plotting	
    #hsc         = (26, 30)    # color scale limits for horizontal plots
    #hsc         = (18, 26)    # good for 100 m
    #hsc         = (14, 19)    # good for 200 m
    #hsc         = (7, 11)     # good for 400 m
    #hsc         = (3.5, 4.4)      # good for 2000 m
    hsc         = (34.4, 34.8)
    d           = 3           # subsampling factor for decent quiver plot [integer]
    sc          = 5           # scale to quiver arrows

# PLOT = 2, 3, 5, 6 ==============================================================

ni          = -15        # latitude, in case PLOT = 2
nj          = -30         # longitude, in case PLOT = 3 
dz          = 10          # delta z [m] for S -- > Z vertical interpolation
Xlim        = (-39, -35)  # x axis limits
#Xlim        = (-23, -10)
Zlim        = (-5000, 0)  # z axis limits
vsc         = (-0.35, 0.35+0.03) # color scale limits for vertical plot
# vsc         = (0, 31) # temp
# vsc         = (34.4, 37.2) # salt
vst         = 0.03        # vel
# vst         = 2         # temp
# vst         = 0.1       # salt

# PLOT = 4 =======================================================================

if PLOT == 4:
	lonplot     = (-41, -32)  # longitude limits for plotting
	latplot     = (-23, -10)  # latitude limits for plotting
	hsc         = (-0.1, 0.1)    # color scale limits for horizontal plots
	d           = 4           # subsampling factor for decent quiver plot [integer]
	sc          = 3           # scale to quiver arrows


# In case you want to loop between time steps, vertical levels
# or different coordinates do create slices, here is the time
# to do it. Below, there is an example on how to do it for batch
# plotting of a sequence of model time steps.  

Tv = []; zetam = [];

if movie == 0: 
    printstr = 0
    execfile('plot_roms.py')
else:
    printstr = 1
    for l in days:
        print " DAY = " + str(l+1)
        execfile('plot_roms.py')
        plt.savefig("figures/surface_vel_day%s.pdf" %(l) )
        plt.clf()
        #tv = transp(x, z, prop)
	#Tv = np.hstack((Tv, tv))
	#zetam = np.hstack((zetam, ZETA.mean()))
	



