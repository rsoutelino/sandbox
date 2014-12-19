#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import seawater.csiro as sw
import scipy.io as sp
from romslab import near2d, near

################################################################################

class FeatureModel(object):
    def __init__(self, filename):
        self.filename = filename
        self.matfile = sp.loadmat(filename)
        self.varlist = self.matfile.keys()   
        for var in self.varlist:
            exec("self.%s = self.matfile.pop('%s')" %(var, var))  
            
def brunt_vaissala(rho, depth):
    """
    Computes Brunt-Vaisalla frequency
    n2 = brunt_vaissala(rho)
        rho: rho profile [1D-array] 
        depth: depth [1D-array]
    """    
    drho = np.gradient(rho)
    dz   = np.gradient(depth)
    g    = 9.8
    rho0 = 1024
    N2   = (g / rho0) * (drho / dz)
    return N2

def burger(N2, H, f, R):
    """
    Computes Burger Number based on the ratio between baroclinic deformation
    radius and curvature radius of a promontory
    USAGE: Bu = burger(N2, H, f, R)
    INPUT:
        N2: brunt vaissalla frequency based on a mean rho profile of the jet
        H:  undisturbed water depth
        f:  coriolis parameter
        R:  radius of curvature of the promontory
    """
    Bu = (N2*H**2) / (f**2 * R**2)   
    return Bu
    
################################################################################     

filename = "/home/rsoutelino/phd/fm/BC-NBUC_FM.mat"
fm = FeatureModel(filename)
depth = np.arange(0, 1520, 20)


# getting the rho profile in the core of the jets at 19 S
lon0, lat0 = -37.411, -19.015 # as result from a ginput
line, col = near2d(fm.lon, fm.lat, lon0, lat0)

temp = fm.temp[:, line, col]
salt = fm.salt[:, line, col]
rho  = sw.dens(salt, temp, depth)
N2   = brunt_vaissala(rho, depth)
H    = depth.max()
f    = sw.cor(lat0)
R    = 100000 # abrolhos bank: ~100km
Bu   = burger(N2, H, f, R)

#plt.figure()
#plt.plot(Bu, -depth)
#plt.grid()
#plt.title("Burger Number Profile")
#plt.show()

fBC = near(depth, 200)[0][0]
BuBC = Bu[:fBC].mean()
BuNBUC = Bu[fBC:].mean()

print "\n The Burger Number for BC is: %s \n" %str(BuBC)
print "\n The Burger Number for NBUC is: %s \n" %str(BuNBUC)
































