#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################
#  Script  plot_roms.py -- Visualization of ROMS fields
#  Rafael Soutelino - rsoutelino@gmail.com
#  Last Modification: Aug, 2010
###################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4

# classes and functions to the computings
from roms_setup import run_setup, get_depths, zlevs, near

def transp(lon, z, v):
	"""
	Slices some 3D field within some lon, lat limits
	Be carefull with the aproximation on distance computing
	"""
	dx = np.diff(lon, axis=1) * 111 * 1000  # valid only for low latitudes!!!
	aux  = dx[:,0]; aux.shape = (np.size(aux), 1)
	dx = np.concatenate( (dx, aux), axis=1) 

	dz = np.diff(z, axis=0)
	aux  = dz[0,:]; aux.shape = (1, np.size(aux))
	dz = np.concatenate( (dz, aux), axis=0)

	transp = np.abs(dx) * np.abs(dz) * v; transp = transp.sum()
	transp = transp / 1e6

	return transp
	

# PYTHON SCRIPT START ######################################################

print ' \n' + '==> ' + '  READING ASCII METADATA FILE ...\n' + ' '
roms = run_setup(expt + '_run.setup')

print ' \n' + '==> ' + '  READING GRID NETCDF FILE ...\n' + ' '
grdfile  = netCDF4.Dataset(roms.rundir + roms.run_name + '_grd.nc')

# assigning some variables from grid file
lon   = grdfile.variables['lon_rho'][:]
lat   = grdfile.variables['lat_rho'][:]
latu  = grdfile.variables['lat_u'][:]
latv  = grdfile.variables['lat_v'][:]
lonu  = grdfile.variables['lon_u'][:]
lonv  = grdfile.variables['lon_v'][:]
h     = grdfile.variables['h'][:]

print ' \n' + '==> ' + '  READING CHOSEN ROMS OUTPUT NETCDF FILE ...\n' + ' '
outfile  = netCDF4.Dataset(roms.rundir + roms.run_name + '_' + filetype + '.nc')

# assigning some variables from ROMS output file
if PLOT == 4:
	ZETA = outfile.variables['zeta'][l,...]
	UBAR = outfile.variables['ubar'][l,...]
	VBAR = outfile.variables['vbar'][l,...]
else:
	T    = outfile.variables['salt'][l,...]
	S    = outfile.variables['salt'][l,...]
	U    = outfile.variables['u'][l,...]
	V    = outfile.variables['v'][l,...]
	print ' \n' + '==> ' + '  GETTING DEPTHS OF S-LEVELS ...\n' + ' '
	km, im, jm = T.shape
	zt   = get_depths(outfile,grdfile,l,'temp')
	zu   = get_depths(outfile,grdfile,l,'u')
	zv   = get_depths(outfile,grdfile,l,'v')

if PLOT == 1 or PLOT == 4:

	# setting up map projection
	m = Basemap(projection='merc',
		llcrnrlat = latplot[0], urcrnrlat = latplot[1], 
		llcrnrlon = lonplot[0], urcrnrlon = lonplot[1],
		 lat_ts=0, resolution='i')
	
	if PLOT == 1:
		print ' \n' + '==> ' + '  INTERPOLATING FROM S --> Z ...\n' + ' '

		u = 0*lonu; v = 0*lonv; t = 0*lon; s = 0*lon

		for a in range (0, im):
			for b in range(0, jm):
				t[a,b] = np.interp(-nk, zt[:, a, b], T[:, a, b] )
				s[a,b] = np.interp(-nk, zt[:, a, b], S[:, a, b] )

		for a in range (0, im):
			for b in range(0, jm-1):
				u[a,b] = np.interp(-nk, zu[:, a, b], U[:, a, b] )

		for a in range (0, im-1):
			for b in range(0, jm):
				v[a,b] = np.interp(-nk, zv[:, a, b], V[:, a, b] )
		
	print ' \n' + '==> ' +\
		'  INTERPOLATING VECTORIAL VARIABLES TO SCALAR RHO POINTS ...\n' + ' '	

	if PLOT == 1:
		u = griddata(lonu.ravel(), latu.ravel(), u.ravel(), lon, lat) 
		v = griddata(lonv.ravel(), latv.ravel(), v.ravel(), lon, lat)

	if PLOT == 4:
		ubar = griddata(lonu.ravel(), latu.ravel(), UBAR.ravel(), lon, lat) 
		vbar = griddata(lonv.ravel(), latv.ravel(), VBAR.ravel(), lon, lat)

	# masking in land values
	if PLOT == 1:
		t = np.ma.masked_where(h < nk, t)
		s = np.ma.masked_where(h < nk, s)
		u = np.ma.masked_where(h < nk, u)
		v = np.ma.masked_where(h < nk, v)

	print ' \n' + '==> ' + '  PLOTTING FIGURE ...\n' + ' '

	lonm, latm = m(lon, lat)
	if PLOT == 1:
		up   = u[::d, ::d]
		vp   = v[::d, ::d]
		prop = t.copy()
	
	if PLOT == 4:
		ZETA = np.ma.masked_where(ZETA > 30, ZETA)
		up   = ubar[::d, ::d]
		vp   = vbar[::d, ::d]
		prop = ZETA.copy()

	lonp = lonm[::d, ::d]
	latp = latm[::d, ::d]

	up   = np.ma.masked_where(up > 1e30, up)
	vp   = np.ma.masked_where(vp > 1e30, vp)
	prop = np.ma.masked_where(prop > 1e30, prop)
	
	fig1 = plt.figure(1, figsize=(6,8), facecolor='w')
	if PLOT == 1:
		#m.pcolormesh(lonm,latm,prop, vmin=hsc[0], vmax=hsc[1],
			#cmap=plt.cm.Spectral_r)
		#plt.colorbar(shrink=0.8)
		m.quiver(lonp, latp, up*sc, vp*sc, scale=10)
		
	if PLOT == 4:
		m.pcolormesh(lonm,latm,prop, vmin=hsc[0], vmax=hsc[1],
			cmap=plt.cm.RdBu_r)
		plt.colorbar(shrink=0.8)
#		m.contour(lonm,latm,prop,20, colors='k')
                m.quiver(lonp, latp, up*sc, vp*sc, scale=10)

	zb = np.ma.masked_where(h >= 100, h)
	m.contourf(lonm, latm, zb, colors=('0.7'), alpha=0.5) 
	text = "20 cm/s";
	ax,ay = m( lonplot[0] + 0.5, latplot[1] - 3.1 ) 
	m.quiver(ax,ay,0.2*sc,0,scale=10,color='k',zorder=10)
	ax,ay = m(lonplot[0] + 0.5, latplot[1] - 3 ); 
	plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
	m.drawcoastlines(zorder=5)
	m.drawcountries(zorder=4)
	m.fillcontinents(color=('0.8'), lake_color='aqua', zorder=3)
	m.drawparallels(np.arange(latplot[0], latplot[1], 1),
		labels=[1, 0, 0, 0], dashes=[1,1000], zorder=6)
	m.drawmeridians(np.arange(lonplot[0], lonplot[1], 2),
		labels=[0, 0, 0, 1], dashes=[1,1000], zorder=7)

	tx, ty = m( lonplot[0] + 0.5, latplot[1] - 1 )
	text = 'Day: ' + str(l+1)
	plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')
	
	if PLOT == 1:
		tx, ty = m( lonplot[0] + 0.5, latplot[1] - 2 )
		text = 'Z: '+ str(nk) +' m' 
		plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')

	if PLOT == 4:
		tx, ty = m( lonplot[0] + 0.5, latplot[1] - 2 )
		text = '$\eta$ x Vbar' 
		plt.text(tx,ty,text,color='k',fontsize=10,fontweight='bold')

	if printstr == 1:
		print "Saving PNG figure ......"
		if PLOT == 1:
			string =  "plt.savefig('figures/"+ expt + filetype +"_vel-salt_" \
					+ str(l+1) +"_"+ str(nk) +"m.png')"
		
		if PLOT == 4:
			string =  "plt.savefig('figures/"+ expt + filetype +"_vbar-zeta_" \
					+ str(l+1) +".png')"

		exec(string)
	else:
		plt.show()

else: 

	zi = dz * np.ceil( zt.min() / dz )
	zi = np.arange(zi, 0+dz, dz)

	print ' \n' + '==> ' + '  INTERPOLATING FROM S --> Z ...\n' + ' '


	if PLOT == 2:
		xsec = latv[:,0]
		fsec = near(xsec,ni)
	   	xsec = lonv[0,:]
	   	z    = np.squeeze( zv[:, fsec, :] )
	   	prop = np.squeeze( V[: ,fsec , :] ) 
  		
	elif PLOT == 3:  
		xsec = lonu[0, :]
		fsec = near(xsec,nj)
	   	xsec = latu[:, 0]
	   	z    = np.squeeze( zu[:, :, fsec] )
	   	prop = np.squeeze( U[: ,: , fsec] ) 

        elif PLOT == 5:
		xsec = latv[:,0]
		fsec = near(xsec,ni)
	   	xsec = lonv[0,:]
	   	z    = np.squeeze( zt[:, fsec, :] )
	   	prop = np.squeeze( T[: ,fsec , :] ) 
  		
	elif PLOT == 6:  
		xsec = lonu[0, :]
		fsec = near(xsec,nj)
	   	xsec = latu[:, 0]
	   	z    = np.squeeze( zt[:, :, fsec] )
	   	prop = np.squeeze( T[: ,: , fsec] ) 
	   	
	xsec.shape = (1, xsec.size)
	x = np.repeat(xsec, roms.klevels, axis=0)
	prop = np.ma.masked_where(prop > 1e10, prop)

	fig1 = plt.figure(1, figsize=(12,6), facecolor='w')
	ax = plt.axes(axisbg=('k'))

        if PLOT == 2 or PLOT == 3:
	    plt.contourf(x,z,prop, np.arange(vsc[0], vsc[1], vst),
		cmap=plt.cm.RdBu)
        else:
            plt.contourf(x,z,prop, np.arange(vsc[0], vsc[1], vst),
		cmap=plt.cm.Spectral_r)
        
	plt.colorbar()
	plt.axis([Xlim[0], Xlim[1], Zlim[0], Zlim[1]])
	if PLOT == 2:	
		plt.xlabel('Longitude')
	else:
		plt.xlabel('Latitude')
	plt.ylabel('Depth')
	plt.title('Cross-section Velocity - Day: '+ str(l+1))

	if printstr == 1:
		print "Saving PNG figure ......"
		string =  "plt.savefig('figures/"+ expt + filetype +"_vel-sec_" \
					+ str(l+1) +".png')"
		exec(string)
	else:
		plt.show()

