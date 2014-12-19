#!/usr/bin/env python
###################################################################################
#  Script  plot_sec_hycom.py -- Visualization of HYCOM fields
#  HYCOM + NCODA Global 1/12 Analysis (expt_60.5) (assimilates data) Nov-2003 to Dec-2006
#  url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_60.5';
#  
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Sep, 2010
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from matplotlib.mlab import find
import datetime as dt
import netCDF4 as nc


print 'Loading U, V, T, S arrays for each time-step and computing mean field:'

U = 0; V = 0; T = 0; S = 0
c = 0

while c <= 17:
	print '    Loading netCDF part: '+ str(c)
	ncfile = nc.Dataset('outputs/hycom_glba0.08_exp60.5_2005_part'+ str(c) +'.nc',
			mode='r', format='NETCDF3_CLASSIC')

	U += ncfile.variables['u'][:]
	V += ncfile.variables['v'][:]
	T += ncfile.variables['temp'][:]
	S += ncfile.variables['salt'][:]
	
	c = c + 1
	del ncfile

ncfile = nc.Dataset('outputs/hycom_glba0.08_exp60.5_2005_part1.nc',
			mode='r', format='NETCDF3_CLASSIC')

lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
z   = ncfile.variables['depth'][:]

U = U / c
V = V / c
T = T / c
S = S / c

U = np.ma.masked_where(U > 1e20, U)
V = np.ma.masked_where(V > 1e20, V)
T = np.ma.masked_where(T > 1e20, T)
S = np.ma.masked_where(S > 1e20, S)

u = np.squeeze(U[0,:,:])
v = np.squeeze(V[0,:,:])
t = np.squeeze(T[0,:,:])
s = np.squeeze(S[0,:,:])

# subsampling:
d = 2

ud   = u[1::d, 1::d]
vd   = v[1::d, 1::d]
lond = lon[1::d, 1::d]
latd = lat[1::d, 1::d]

plt.figure()
plt.pcolormesh(lon, lat, t, cmap=plt.cm.Spectral_r); plt.colorbar()
plt.quiver(lond, latd, ud, vd, scale=5)
