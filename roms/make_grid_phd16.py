#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
#
#  Build a ROMS grid file
#
#  Further Information:  
#  http://www.brest.ird.fr/Roms_tools/
#  
#  This file is part of ROMSTOOLS
#
#  ROMSTOOLS is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2002-2006 by Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#
#  Contributions of P. Marchesiello (IRD) and X. Capet (UCLA)
#
#  Updated    Aug-2006 by Pierrick Penven
#  Updated    24-Oct-2006 by Pierrick Penven (mask correction)
#
#  Translated to Python by Rafael Soutelino: rsoutelino@gmail.com 
#  Last Modification: Aug, 2010
################################################################


print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 
# IMPORTING MODULES #################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import delaunay
from matplotlib.mlab import griddata, find
from mpl_toolkits.basemap import Basemap
from scipy.io import loadmat
import datetime as dt
import netCDF4

# classes and functions to the computings
from roms_setup import run_setup, rho2uvp, get_metrics, spheric_dist
from roms_setup import get_angle, add_topo, process_mask, uvp_mask, smoothgrid
from roms_setup import rotfilter, rfact, hanning_smoother
from roms_setup import hanning_smoother_coef2d, FX, FY 

# SCRIPT START ######################################################

# Basic Settings:

filenamestr = '_grd.nc'
filetypestr = 'ROMS Grid file'

# READING PREVIOUSLY BUILT RELEVANT FILES: ###########################
# metadata ascii file
# OA-created netcdf initial T, S file 
# grid netcdf file

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
run = run_setup('../phd16_run.setup')

dl = run.resolution
lonr  = np.arange(run.lonmin, run.lonmax + dl, dl)

i = 0; 
latr = np.array([run.latmin])
while latr[i] <= run.latmax:
    i    = i + 1
    tmp  = latr[i-1] + dl * np.cos( latr[i-1]*np.pi/180 )
    latr = np.hstack([latr, tmp])



Lonr, Latr        =  np.meshgrid(lonr, latr)
Lonu, Lonv, Lonp  = rho2uvp(Lonr)
Latu, Latv, Latp  = rho2uvp(Latr)

M, L = Latp.shape

print ' \n' + '==> ' + '  COMPUTING METRICS  ...\n' + ' '


print ' \n' + '==> ' + '  LLm = ' + np.str(L-1) + ' ...\n' + ' '
print ' \n' + '==> ' + '  MMm = ' + np.str(M-1) + ' ...\n' + ' '

# !!!!!!!!!!!!!!!!!!!!!
### CODE SOMETHING HERE TO WRITE THIS INFORMATION IN THE METADATA FILE
# !!!!!!!!!!!!!!!!!!!!!

pm, pn, dndx, dmde = get_metrics(Latu, Lonu, Latv, Lonv)
xr = 0*pm
yr = xr

for i in np.arange(0, L):
    xr[:, i+1] = xr[:, i] + 2 / ( pm[:, i+1] + pm[:, i] )

for j in np.arange(0, M):
    yr[j+1, :] = yr[j, :] + 2 / ( pn[j+1, :] + pn[j, :] )

xu, xv, xp = rho2uvp(xr)
yu, yv, yp = rho2uvp(yr)

dx    = 1 / pm
dy    = 1 / pn
dxmax = np.max( dx/1000 )
dxmin = np.min( dx/1000 )
dymax = np.max( dy/1000 )
dymin = np.min( dy/1000 )

angle = get_angle(Latu, Lonu)

f0 = 4 * np.pi * np.sin( np.pi * Latr/180 ) / ( 24*3600 )


print ' \n' + '==> ' + '  ADDING TOPOGRAPHY ...\n' + ' '

h = add_topo(Lonr, Latr, pm, pn, run.topo_filename)
hraw = h.copy()
h[ np.where(h > run.hmax) ] = run.hmax

print ' \n' + '==> ' + '  COMPUTING THE MASK ...\n' + ' '

maskr = h*0
maskr[ np.where(h > 0) ] = 1 
maskr = process_mask(maskr)
[masku, maskv, maskp] = uvp_mask(maskr) 

print ' \n' + '==> ' + '  FILTERING THE TOPOGRAPHY ...\n' + ' '

h = smoothgrid(h, maskr, run.hmin, run.hmaxc,
             run.slope, run.npass, run.nfinal)


####################################################################
####################################################################

print ' \n' + '==> ' + '  WRITING NETCDF GRID FILE ...\n' + ' '

now = dt.datetime.now()
Lp = L + 1
Mp = M + 1

if run.spherical == 1:
	spherical = 'T'
else:
	spherical = 'F'

ncfile = netCDF4.Dataset(run.rundir + run.run_name + filenamestr, mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_psi', size=L)
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_psi', size=M)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
#ncfile.createDimension('bath', size=None)
ncfile.createDimension('one', size=1)
ncfile.createDimension('two', size=2)
ncfile.createDimension('four', size=4)


# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', filetypestr)
setattr(ncfile, 'title', run.topo_info)
setattr(ncfile, 'history', str(now))


# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
setattr(ncfile.variables['spherical'], 'option_T', 'spherical')
setattr(ncfile.variables['spherical'], 'option_F', 'cartesian')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('xl', 'd', dimensions=('one'))
setattr(ncfile.variables['xl'], 'long_name', 'domain length in XI-direction')
setattr(ncfile.variables['xl'], 'units', 'meter')
ncfile.variables['xl'][:]  = xr.max()

# ---------------------------------------------------------------------------
ncfile.createVariable('el', 'd', dimensions=('one'))
setattr(ncfile.variables['el'], 'long_name', 'domain length in ETA-direction')
setattr(ncfile.variables['el'], 'units', 'meter')
ncfile.variables['el'][:]  = yr.max()

# ---------------------------------------------------------------------------
ncfile.createVariable('hraw', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['hraw'], 'long_name', 'Working bathymetry at RHO-points')
setattr(ncfile.variables['hraw'], 'units', 'meter')
setattr(ncfile.variables['hraw'], 'coordinates', 'lon_rho lat_rho bath')
ncfile.variables['hraw'][:]  = hraw

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'Final bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:]  = h

# ---------------------------------------------------------------------------
ncfile.createVariable('f', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['f'], 'long_name', 'Coriolis parameter at RHO-points')
setattr(ncfile.variables['f'], 'units', 'second-1')
setattr(ncfile.variables['f'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['f'][:]  = f0

# ---------------------------------------------------------------------------
ncfile.createVariable('pm', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pm'], 'long_name', 'Curvilinear coordinate metric in XI')
setattr(ncfile.variables['pm'], 'units', 'meter-1')
setattr(ncfile.variables['pm'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pm'][:]  = pm

# ---------------------------------------------------------------------------
ncfile.createVariable('pn', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pn'], 'long_name', 'Curvilinear coordinate metric in ETA')
setattr(ncfile.variables['pn'], 'units', 'meter-1')
setattr(ncfile.variables['pn'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pn'][:]  = pn

# ---------------------------------------------------------------------------
ncfile.createVariable('dndx', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dndx'], 'long_name', 
	'XI derivative of inverse metric factor pn')
setattr(ncfile.variables['dndx'], 'units', 'meter')
setattr(ncfile.variables['dndx'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dndx'][:]  = dndx

# ---------------------------------------------------------------------------
ncfile.createVariable('dmde', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dmde'], 'long_name', 
	'ETA derivative of inverse metric factor pm')
setattr(ncfile.variables['dmde'], 'units', 'meter')
setattr(ncfile.variables['dmde'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dmde'][:]  = dmde

# ---------------------------------------------------------------------------
ncfile.createVariable('x_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['x_rho'], 'long_name', 'x location of RHO-points')
setattr(ncfile.variables['x_rho'], 'units', 'meter')
ncfile.variables['x_rho'][:]  = xr

# ---------------------------------------------------------------------------
ncfile.createVariable('y_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['y_rho'], 'long_name', 'y location of RHO-points')
setattr(ncfile.variables['y_rho'], 'units', 'meter')
ncfile.variables['y_rho'][:]  = yr

# ---------------------------------------------------------------------------
ncfile.createVariable('x_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['x_psi'], 'long_name', 'x location of PSI-points')
setattr(ncfile.variables['x_psi'], 'units', 'meter')
ncfile.variables['x_psi'][:]  = xp

# ---------------------------------------------------------------------------
ncfile.createVariable('y_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['y_psi'], 'long_name', 'y location of PSI-points')
setattr(ncfile.variables['y_psi'], 'units', 'meter')
ncfile.variables['y_psi'][:]  = yp

# ---------------------------------------------------------------------------
ncfile.createVariable('x_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['x_u'], 'long_name', 'x location of U-points')
setattr(ncfile.variables['x_u'], 'units', 'meter')
ncfile.variables['x_u'][:]  = xu

# ---------------------------------------------------------------------------
ncfile.createVariable('y_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['y_u'], 'long_name', 'y location of U-points')
setattr(ncfile.variables['y_u'], 'units', 'meter')
ncfile.variables['y_u'][:]  = yu

# ---------------------------------------------------------------------------
ncfile.createVariable('x_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['x_v'], 'long_name', 'x location of V-points')
setattr(ncfile.variables['x_v'], 'units', 'meter')
ncfile.variables['x_v'][:]  = xv

# ---------------------------------------------------------------------------
ncfile.createVariable('y_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['y_v'], 'long_name', 'y location of V-points')
setattr(ncfile.variables['y_v'], 'units', 'meter')
ncfile.variables['y_v'][:]  = yv

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree east')
ncfile.variables['lon_rho'][:]  = Lonr

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree north')
ncfile.variables['lat_rho'][:]  = Latr

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lon_psi'], 'long_name', 'longitude of PSI-points')
setattr(ncfile.variables['lon_psi'], 'units', 'degree east')
ncfile.variables['lon_psi'][:]  = Lonp

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lat_psi'], 'long_name', 'latitude of PSI-points')
setattr(ncfile.variables['lat_psi'], 'units', 'degree north')
ncfile.variables['lat_psi'][:]  = Latp

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree east')
ncfile.variables['lon_u'][:]  = Lonu

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree north')
ncfile.variables['lat_u'][:]  = Latu

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree east')
ncfile.variables['lon_v'][:]  = Lonv

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree north')
ncfile.variables['lat_v'][:]  = Latv

# ---------------------------------------------------------------------------
ncfile.createVariable('angle', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['angle'], 'long_name', 'angle between XI-axis and EAST')
setattr(ncfile.variables['angle'], 'units', 'radians')
setattr(ncfile.variables['angle'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['angle'][:]  = angle

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['mask_rho'], 'long_name', 'mask on RHO-points')
setattr(ncfile.variables['mask_rho'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_rho'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_rho'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_rho'][:]  = maskr

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['mask_u'], 'long_name', 'mask on U-points')
setattr(ncfile.variables['mask_u'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_u'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_u'], 'coordinates', 'lon_u lat_u')
ncfile.variables['mask_u'][:]  = masku

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['mask_v'], 'long_name', 'mask on V-points')
setattr(ncfile.variables['mask_v'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_v'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_v'], 'coordinates', 'lon_v lat_v')
ncfile.variables['mask_v'][:]  = maskv

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['mask_psi'], 'long_name', 'mask on PSI-points')
setattr(ncfile.variables['mask_psi'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_psi'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_psi'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_psi'][:]  = maskp

ncfile.sync()

print ' \n' + '==> ' + '  ############################################  ...\n' + ' '
print ' \n' + '==> ' + '        GRID FILE SUCCESSFULLY CREATED          ...\n' + ' '
print ' \n' + '==> ' + '  ############################################  ...\n' + ' '















































































