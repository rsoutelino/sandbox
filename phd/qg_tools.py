#!/usr/local/epd-7.2-1-rh5-x86_64/bin/python
# -*- coding:utf-8 -*-

import numpy as np
from matplotlib.mlab import griddata
import seawater.csiro as sw
import scipy.io as sp
import romslab

################################################################################

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

def compute_burger(N2, H, f, R):
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
    
def compute_dqds(N2, f, v, depth):
    """
    Computes mean potential vorticity cross-jet gradient (dq/ds)
    USAGE: dqds = dqds(N2, f, v)
    INPUT:
        N2: brunt vaissalla frequency based on a mean rho profile of the jet
        f:  coriolis parameter
        v:  along jet velocity
    """
    part1 = f**2 / N2
    part2 = np.gradient(v) / np.gradient(depth)
    part3 = part1 * part2
    dqds  = np.gradient(part3) / np.gradient(depth)
    return dqds

