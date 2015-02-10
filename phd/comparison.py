# -*- coding: utf-8 -*-
######################################################
## plotando os dados do ADCP processados pelo CODAS ##
##         Rafael Soutelino - set/2007              ##
######################################################
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
from mpl_toolkits.basemap import Basemap
from matplotlib.mlab import griddata
import romslab
import netCDF4 as nc

# defining functions

def resh(x):
    """
    Reshapes 2D arrays to 1D
    """
    a = x.shape[0]
    b = x.shape[1]
    return x.reshape(a*b, 1), a, b

def near(x,x0):
    import pylab as pl 
    dx = x - x0
    dx = abs(dx)
    fn = pl.find(dx == dx.min())
    return fn

def mask_nans(array):
    f = np.where(np.isnan(array)==1)
    array[f] = -99999
    array = np.ma.masked_where(array==-99999, array)
    return array
    
#################################################
# preparing to plot figures
lims = [-41, -33.7, -22.5, -12]
print "Basemap"
m = Basemap(projection='cyl', llcrnrlat=-22.5, urcrnrlat=-12.0, llcrnrlon=-41.0, 
            urcrnrlon=-33.7, lat_ts=0, resolution='i')
s = 2; # escala para vetores
ct = 200; # ct2 = 20;
pp = 1000;
p = 50
etopo = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/etopo2_leste.mat')
xb = etopo.pop('xb')
yb = etopo.pop('yb')
zb = etopo.pop('zb')
xb,yb = m(xb,yb)
# lon,lat = m(lon,lat)
# masking topography to plot continental shelf mask\
zb = pl.ma.masked_where(zb <= -pp, zb)


#################################################
fig1 = plt.figure(1,figsize=(15,8),facecolor='w')

###################################################################################

oe2 = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/psiobs_OEII.mat')
psiob = oe2.pop('psiob')
xgi = oe2.pop('xgi'); ygi = oe2.pop('ygi');
uo = oe2.pop('uo'); vo = oe2.pop('vo');  

s = 2;

# masking 
f = pl.find(pl.isnan(psiob) == 1); 
a, b = psiob.shape; 
psiob = pl.ravel(psiob); uo = pl.ravel(uo); vo = pl.ravel(vo);

psiob[f] = 0; uo[f] = 0; vo[f] = 0;
psiob.shape = (a,b); uo.shape = (a,b); vo.shape = (a,b)
psiob = pl.ma.masked_where(psiob == 0, psiob)
uo = pl.ma.masked_where(uo == 0, uo)
vo = pl.ma.masked_where(vo == 0, vo)

g = 1

p1 = plt.subplot(131)
#con = m.contourf(xgi,ygi,psigo,scale,linewidth=0,alpha = 0.5); 
m.quiver(xgi[::g,::g], ygi[::g,::g], uo[::g,::g]*s, vo[::g,::g]*s, scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
  #[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
  #set(h,'facecolor',[.7 .7 .7]);
  #[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
  #set(h,'color',[.7 .7 .7]);
m.drawparallels(np.arange(-22,-12,2),labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-38.5,-18.5); plt.text(ax,ay,"AB",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VTR",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p1.text(ax, ay, "Mar 2005", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p1.text(ax, ay, "OEII ADCP-derived", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p1.text(ax, ay, "a", ha="center", va="center", size=12,bbox=bbox_props)

######################################################################

# get roms velocity for day 30 of CONTROL run
outname = "/home/rsoutelino/myroms/phd_run/phd15_avg.nc"
roms = romslab.RomsHis(outname)
lims = [-41, -34.5, -20, -12]
U, V = roms.u[30, ...], roms.v[30,...]
zlevs = romslab.get_depths(roms.ncfile, 0, 'temp')
lon, lat, U = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], U, lims)
lon, lat, V = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], V, lims)
lon, lat, zlevs = romslab.subset(roms.lon_rho[:], roms.lat_rho[:], zlevs, lims)
u, v = romslab.get_hslice(50, zlevs, U), romslab.get_hslice(50, zlevs, V)
u, v = mask_nans(u), mask_nans(v)
u, v = np.ma.masked_where(u > 10, u), np.ma.masked_where(v > 10, v)

xgi = xgi[:-10, :]; ygi = ygi[:-10, :];

u = griddata(lon.ravel(), lat.ravel(), u.ravel(), xgi, ygi)
v = griddata(lon.ravel(), lat.ravel(), v.ravel(), xgi, ygi)

p2 = plt.subplot(132)
#m.contourf(xgi,ygi,psiob,scale,linewidth=0,alpha = 0.5);
m.quiver(xgi[::g,::g], ygi[::g,::g], u[::g,::g]*s, v[::g,::g]*s, scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[0, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-38.5,-18.5); plt.text(ax,ay,"AB",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VTR",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p2.text(ax, ay, "Day 30", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p2.text(ax, ay, "ROMS CONTROL Run", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p2.text(ax, ay, "b", ha="center", va="center", size=12,bbox=bbox_props)

######################################################################

# OCCAM

s = 2; # escala para vetores
ct = 200; # ct2 = 20;
pp = 1000;
p = 50
xgi,ygi = m(xgi, ygi)
occam = nc.Dataset('/home/rsoutelino/phd/occam/data/y2003.h5m1subvol')
lon = occam.variables['LONGITUDE_U'][:] - 360
lat = occam.variables['LATITUDE_U'][:]
lon, lat = np.meshgrid(lon, lat)
dataset = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/occam_feb2003.mat')
u   = dataset.pop('u')  ; v   = dataset.pop('v')  ;

u = griddata(lon.ravel(), lat.ravel(), u.ravel(), xgi, ygi)
v = griddata(lon.ravel(), lat.ravel(), v.ravel(), xgi, ygi)

p3 = plt.subplot(133)

m.quiver(xgi, ygi, u*s, v*s, scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[0, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-38.5,-18.5); plt.text(ax,ay,"AB",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VTR",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p3.text(ax, ay, "Feb 2003", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p3.text(ax, ay, "OCCAM", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p3.text(ax, ay, "c", ha="center", va="center", size=12,bbox=bbox_props)


#############################################################################################
fig1.tight_layout(pad=3)
plt.show()
plt.savefig('figures/comparison.pdf')


  


