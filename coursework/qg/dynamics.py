#!/usr/local/epd-7.2-1-rh5-x86_64/bin/python
# -*- coding:utf-8 -*-
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import seawater.csiro as sw
import scipy.io as sp
import romslab
from qg_tools import *
from qg_plot import *

################################################################################

class FeatureModel(object):
    def __init__(self, filename):
        self.filename = filename
        self.matfile = sp.loadmat(filename)
        self.varlist = self.matfile.keys() 
        self.depth = np.arange(0, 5020, 20)
        for var in self.varlist:
            exec("self.%s = self.matfile.pop('%s')" %(var, var))  

def mask_nans(array):
    f = np.where(np.isnan(array)==1)
    array[f] = -99999
    array = np.ma.masked_where(array==-99999, array)
    return array
            
################################################################################     

mainpath = "/home/rsoutelino/"

filename = mainpath + "phd/fm/BC-NBUC_FM.mat"
fm = FeatureModel(filename)


# BURGER NUMBER ################

# getting the rho profile in the core of the jets at 19 S
lon0, lat0 = -37.411, -19.015 # as result from a ginput
line, col = romslab.near2d(fm.lon, fm.lat, lon0, lat0)

temp = fm.temp[:, line, col]
salt = fm.salt[:, line, col]
depth = fm.depth
rho  = sw.dens(salt, temp, depth)
N2   = brunt_vaissala(rho, depth)
H    = depth.max()
f    = sw.cor(lat0)
R    = 100000 # abrolhos bank: ~100km (the most wide one)
Bu   = compute_burger(N2, H, f, R)

fBC = romslab.near(depth, 300)[0][0]
BuBC = Bu[:fBC].mean()
BuNBUC = Bu[fBC:].mean()

# ROSSBY NUMBER ################

fNBUC = romslab.near(depth, 400)[0][0]
f0 = sw.cor(fm.lat.mean())
L = 150000 
Ubc   = np.sqrt(fm.u[0,...].max()**2 + fm.v[0,...].max()**2)
Unbuc = np.sqrt(fm.u[fNBUC,...].max()**2 + fm.v[fNBUC,...].max()**2)

RoBC = np.abs(Ubc / (f0*L))
RoNBUC = np.abs(Unbuc / (f0*L))

print "\n The Burger Number for BC-FM is: %s " %str(BuBC)
print "\n The Burger Number for NBUC-FM is: %s " %str(BuNBUC)
print "\n The Rossby Number for BC-FM is: %s " %str(RoBC)
print "\n The Rossby Number for NBUC-FM is: %s " %str(RoNBUC)

# dq/ds ##########################

# getting the sectional distribution at 19 S, from ginput
line, col1 = romslab.near2d(fm.lon, fm.lat, -38.0, -18.5)
line, col2 = romslab.near2d(fm.lon, fm.lat, -36.45, -18.5)
v = fm.v[:, 0, col1:col2].mean(axis=1)
dvdz = np.gradient(v) / np.gradient(depth)
dqdn = compute_dqds(N2, f, v, depth)


ilson = sp.loadmat('ilson.mat')
v     = ilson['v']
dvdz  = ilson['dvdz']
dqdn  = ilson['dqdn']
depth = ilson['depth']

figname = "figures/instability_conditions.pdf"
plot_instability_conditions(v*100, dvdz, dqdn, depth, figname)


figname = "figures/instab_2.pdf"
L     = ilson['L'].ravel()
sigma = ilson['sigma'].ravel()
cp    = ilson['cp'].ravel()
gv    = ilson['gv'].ravel()

plt.figure(facecolor='w')
p1 = plt.subplot(2,1,1)
plt.plot(L, sigma, 'r', linewidth=2)
plt.ylabel('$\sigma$ [days$^{-1}$]')
# plt.xlabel('Wavelength [km]')
plt.grid()
p1.set_axis_bgcolor('0.95')

p2 = plt.subplot(2,1,2)
p = plt.plot(L, cp, 'b', L, gv, 'g--', linewidth=2)
plt.ylabel('Phase speed and group\n velocity [m s$^{-1}$]')
plt.xlabel('Wavelength [km]')
plt.grid()
p2.set_axis_bgcolor('0.95')
p2.set_ylim([-0.005, 0.11])
plt.show()
plt.savefig(figname)
stop

# run-averaged ROMS fields ########

gridname = mainpath + "myroms/phd_run/phd15_grd.nc"
outname = mainpath +  "myroms/phd_run/azores_run/phd15_avg-1.nc"
grid = romslab.RomsGrid(gridname)
roms = romslab.RomsHis(outname)

# limits for the Abrolhos, Royal-Charlotte  and Ilheus Eddies
lims = -40, -36, -19.5, -14.0
# lims = -40, -36, -18, -14.0
M = Basemap(projection='cyl', llcrnrlon=lims[0], urcrnrlon=lims[1],
                llcrnrlat=lims[2], urcrnrlat=lims[3], lat_ts=0, resolution='i')
depths = 0, 400
print "\nComputing run-averaged fields =======================\n\n"
lon, lat, U = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], roms.u, lims)
lon, lat, V = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], roms.v, lims)
lon, lat, T = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], roms.temp, lims)
lon, lat, h     = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], roms.h, lims)
UbarS, VbarS, TbarS = U.mean(axis=0), V.mean(axis=0), T.mean(axis=0)

zlevs = romslab.get_depths(roms.ncfile, 30, 'temp')
lon, lat, zlevs = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], zlevs, lims)
Tbar, Ubar, Vbar = [],[],[]

for k in depths:
    aux = romslab.get_hslice(k, zlevs, TbarS)
    aux = mask_nans(aux)
    Tbar.append( aux )
    aux = romslab.get_hslice(k, zlevs, UbarS)
    aux = mask_nans(aux)
    Ubar.append( aux )
    aux = romslab.get_hslice(k, zlevs, VbarS)
    aux = mask_nans(aux)
    Vbar.append( aux )

figname ="figures/average_fields.pdf"
figname ="figures/average_fields.svg"
plot_average_fields(lon, lat, Ubar, Vbar, VbarS, h, zlevs, depths, lims, M, figname)

# RUN-AVERAGE based Rossby Number     ##########################################
Ro = []
Ro.append( np.abs( ( np.sqrt(Ubar[0]**2 + Vbar[0]**2) ) / (f0*L) ) )
Ro.append( np.abs( ( np.sqrt(Ubar[1]**2 + Vbar[1]**2) ) / (f0*L) ) )
# plot_rossby(lon, lat, M, lims, h, Ro, "figures/rossby_map.pdf")

# BETA PLANE APPROXIMATION #####################################################

df = np.abs( sw.cor( lat.min() ) - sw.cor( lat.max() ) )
dy = np.abs( lat.min() - lat.max() ) * 111000
beta = df / dy
print "\nMean f0 for the region is: %s" %f0
print "BETA for the region is: %s" %beta
print "Y scale is 10^5"
print "f0 >> by"

# computing Instantaneuos Baroclinic Conversions ###############################

fm_bar = False

if fm_bar:
    Ubar = [ romslab.get_hslice(depths[0], zlevs, U[0,...]),
             romslab.get_hslice(depths[1], zlevs, U[0,...]) ]
    Vbar = [ romslab.get_hslice(depths[0], zlevs, V[0,...]),
             romslab.get_hslice(depths[1], zlevs, V[0,...]) ]
    Tbar = [ romslab.get_hslice(depths[0], zlevs, T[0,...]),
             romslab.get_hslice(depths[1], zlevs, T[0,...]) ]
    Tbar = [ mask_nans(array) for array in Tbar ]
    Ubar = [ mask_nans(array) for array in Ubar ]
    Vbar = [ mask_nans(array) for array in Vbar ]


# initializing velocity and temperature perturbations arrays
Unot, Vnot, Tnot = [], [], []
for k in range(len(depths)):
    Unot.append( [ romslab.get_hslice(depths[k], zlevs, U[0,...]) - Ubar[k] ] )
    Vnot.append( [ romslab.get_hslice(depths[k], zlevs, V[0,...]) - Vbar[k] ] )
    Tnot.append( [ romslab.get_hslice(depths[k], zlevs, T[0,...]) - Tbar[k] ] )

# computing all velocity and temperature perturbation arrays
for t in range(1, 35):
    print "Computing perturbations for day %s" %t
    for k in range(len(depths)):
        Unot[k].append( romslab.get_hslice(depths[k], zlevs, U[t,...]) - Ubar[k] )
        Vnot[k].append( romslab.get_hslice(depths[k], zlevs, V[t,...]) - Vbar[k] )
        Tnot[k].append( romslab.get_hslice(depths[k], zlevs, T[t,...]) - Tbar[k] )


# baroclinic conversions #######################################################
BC = [[],[]]
g = 9.8

for k in range(len(depths)):
    theta_z = Tbar[k].mean()
    alpha = sw.alpha( 35, theta_z, depths[k] )
    gradT = np.gradient(Tbar[k]) 
    dTdx = gradT[1] / (0.05 * 111000)
    dTdy = gradT[0] / (0.05 * 111000)
    for t in range(0, 35):
        print "Computing mean baroclinc conversions for day %s" %t
        bc =  ( (g*alpha)/theta_z ) * ( Unot[k][t]*Tnot[k][t]*dTdx + Vnot[k][t]*Tnot[k][t]*dTdy )
        bc = bc.mean()
        BC[k].append(bc)

 
# plot_energy_conversion(BC, "figures/energy_conversion.pdf")
# plt.show()

