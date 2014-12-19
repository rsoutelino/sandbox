# -*- coding: utf-8 -*-
######################################################
## define title!!!!!!!!!!!  date, authorship        ##
######################################################
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sp
from mpl_toolkits.basemap import Basemap
import seawater.csiro as sw

######################################################################
#    FUNCTIONS    ####################################################
######################################################################
print "\n\nDefining Functions %s\n\n" % ("."*50)

def near(x,x0):
    import pylab as pl 
    dx = x - x0
    dx = abs(dx)
    fn = pl.find(dx == dx.min()) 
    return fn

def apply_mask(x):
    """
    Gets arrays with NaN from MAT files and applies python masked_where
    """
    f = pl.find(pl.isnan(x) == 1)
    l1, l2 = x.shape 
    x = pl.ravel(x)
    x[f] = 0
    x.shape = (l1,l2)
    x = pl.ma.masked_where(x == 0, x)
    return x

def subplotfig(sbplt, d, s, name, area, panel, title):
    """
    sbplt : subplot number [integer]
    d     : array sub-sampling factor [integer]
    siz   : figure size [tuple]
    name  : dataset name [string]
    area  : projection sub-area [integer]
    panel : a, b, or c [string] ex: "(c)"
    title : a title to the plot ex: AVISO
    """
    plt.subplot(sbplt)
    if area == 1:
        m = m1
        etopo1['mlon'], etopo1['mlat'] = m( etopo1['lon'],
            etopo1['lat'] )
    else:
        m = m2
        etopo1['mlon'], etopo1['mlat'] = m( etopo1['lon'],
            etopo1['lat'] )
    exec( "m.quiver( %s['mlon'][::d,::d], %s['mlat'][::d,::d], %s['u'][::d,::d]*s, %s['v'][::d,::d]*s, scale=10 )" %(name, name, name, name) )
    if name == 'oe2' or name == 'proab':
        m.contourf( etopo1['mlon'], etopo1['mlat'], etopo1['mh'],
        colors=('w'), linewidth=0, alpha=0.7 )
    m.contour( etopo1['mlon'], etopo1['mlat'], etopo1['h']*-1,
        (1000,200), colors='k' )
    m.drawcoastlines(zorder=5)
    m.drawcountries(zorder=4)
    m.fillcontinents(color='0.0',lake_color='aqua',zorder=3)
    m.drawparallels(np.arange(-25,-12,1),labels=[1, 1, 0, 0],
        dashes=[1,1000],zorder=6)
    if panel == '(a)':
        lab = (0, 0, 1, 1)
    elif panel == '(b)':
        lab = (0, 0, 0, 1)
    else:
        lab = (0, 0, 0, 1)
    m.drawmeridians(np.arange(-44,-32,1),labels=lab,dashes=[1,1000],
        zorder=7)
    if area == 1:
        ax,ay = m(-40.7,-17)
        plt.text(ax, ay, title, color='w', fontsize=10,
            fontweight='bold')
        ax,ay = m(-40.7,-17.5)
        m.quiver(ax, ay, 0.5*s, 0, scale=10, color='w',zorder=10)
        ax,ay = m(-40.7,-17.7)
        plt.text(ax, ay, "50 cm/s", color='w', fontsize=10,
            fontweight='bold')
        ax,ay = m(-40.7,-18.0)
        m.quiver(ax, ay, 0.2*s, 0, scale=10, color='w', zorder=10)
        ax,ay = m(-40.7,-18.2)
        plt.text(ax, ay, "20 cm/s", color='w', fontsize=10, fontweight='bold')
        bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
        ax,ay = m(-40.7,-19)
        plt.text(ax, ay, panel, ha="left", va="center", size=12, bbox=bbox_props)
    else:
        ax,ay = m(-40.7,-17)
        plt.text(ax, ay, title, color='w', fontsize=10,
            fontweight='bold')
        ax,ay = m(-40.7,-17.5)
        m.quiver(ax, ay, 0.5*s, 0, scale=10, color='w',zorder=10)
        ax,ay = m(-40.7,-17.9)
        plt.text(ax, ay, "50 cm/s", color='w', fontsize=10,
            fontweight='bold')
        ax,ay = m(-40.7,-18.5)
        m.quiver(ax, ay, 0.2*s, 0, scale=10, color='w', zorder=10)
        ax,ay = m(-40.7,-18.9)
        plt.text(ax, ay, "20 cm/s", color='w', fontsize=10, fontweight='bold')
        bbox_props = dict(boxstyle="round", fc="w", ec="0.5")
        ax,ay = m(-40.7,-20)
        plt.text(ax, ay, panel, ha="left", va="center", size=12, bbox=bbox_props)
    return name*2 # just to use "return"



######################################################################
#    CONFIGURING    ##################################################
######################################################################
print "\n\nConfiguring Script %s\n\n" % ("."*50)

pathname1 = "/home/rsoutelino/rsoutelino/prod/grl_msc/proc/"
pathname1 = "/media/RAFA500G/rsoutelino/prod/grl_msc/proc/"
pathname2 = ""

filename1 = "etopo2_leste.mat"
filename2 = "psiobs_OEII.mat"
filename3 = "psiobs_proab_10m.mat"

m1 = Basemap(projection='merc',llcrnrlat=-19.5,urcrnrlat=-16.5,
    llcrnrlon=-41.0,urcrnrlon=-35.5,lat_ts=0,resolution='i')
m2 = Basemap(projection='merc',llcrnrlat=-22,urcrnrlat=-16.5,
    llcrnrlon=-41.0,urcrnrlon=-35.5,lat_ts=0,resolution='i')


######################################################################
#    LOADING STUFF    ################################################
######################################################################
print "\n\nLoading Data %s\n\n" % ("."*50)

# Bathymetry =========================================================
pp = 1000
etopo1  = {'lon':[],'lat':[],'mlon':[],'mlat':[],'h':[], 'mh':[]}
dataset = sp.loadmat( "%s%s" % (pathname1, filename1) )
etopo1['lon'] = dataset.pop('xb')
etopo1['lat'] = dataset.pop('yb')
etopo1['h']   = dataset.pop('zb')
etopo1['mh'] = pl.ma.masked_where(etopo1['h'] <= -pp, etopo1['h'])

# LESTE 2 DATA =======================================================
oe2  = {'lon':[],'lat':[],'mlon':[],'mlat':[],'psi':[],'u':[],'v':[]}
dataset = sp.loadmat( "%s%s" % (pathname1, filename2) )
oe2['psi'] = dataset.pop('psiob')
oe2['lon'] = dataset.pop('xgi')
oe2['lat'] = dataset.pop('ygi')
oe2['mlon'], oe2['mlat'] = m1( oe2['lon'], oe2['lat'] )
oe2['u']   = dataset.pop('uo')
oe2['v']   = dataset.pop('vo')
# masking 
oe2['psi'] = apply_mask(oe2['psi'])

# PROABROLHOS DATA ===================================================
proab  = {'lon':[],'lat':[],'mlon':[],'mlat':[],'psi':[],'u':[],'v':[]}
dataset = sp.loadmat( "%s%s" % (pathname2, filename3) )
proab['psi'] = dataset.pop('psiob')
proab['lon'] = dataset.pop('xgi')
proab['lat'] = dataset.pop('ygi')
proab['mlon'], proab['mlat'] = m2( proab['lon'], proab['lat'] )
proab['u']   = dataset.pop('uo')
proab['v']   = dataset.pop('vo')
# masking 
proab['psi'] = apply_mask(proab['psi'])

# AVISO AREA 1 =======================================================
# Hardcoded for pathnames and filenames due to the nature of the 
# MAT files provided
aviso1 = {'lon':[],'lat':[],'mlon':[],'mlat':[],'ssh':[],'u':[], 'v':[]}
lat = sp.loadmat('aviso_files/ay1.mat'); lat = lat.pop('ay1')
lon = sp.loadmat('aviso_files/ax1.mat'); lon = lon.pop('ax1')
lon, lat = np.meshgrid(lon, lat)
aviso1['lat'] = lat
aviso1['lon'] = lon
aviso1['mlon'], aviso1['mlat'] = m1(aviso1['lon'], aviso1['lat'])
ssh = sp.loadmat('aviso_files/assh1.mat'); aviso1['ssh'] = ssh.pop('assh1')
u   = sp.loadmat('aviso_files/au1.mat'); aviso1['u'] = u.pop('au1')
v   = sp.loadmat('aviso_files/av1.mat'); aviso1['v'] = v.pop('av1')
# masking 
aviso1['u']   = apply_mask(aviso1['u'])
aviso1['v']   = apply_mask(aviso1['v'])
aviso1['ssh'] = apply_mask(aviso1['ssh'])
# scaling
aviso1['u'] = aviso1['u']/100; aviso1['v'] = aviso1['v']/100;

# AVISO AREA 2 =======================================================
# Hardcoded for pathnames and filenames due to the nature of the 
# MAT files provided
aviso2 = {'lon':[],'lat':[],'mlon':[],'mlat':[],'ssh':[],'u':[], 'v':[]}
lat = sp.loadmat('aviso_files/ay2.mat'); lat = lat.pop('ay2')
lon = sp.loadmat('aviso_files/ax2.mat'); lon = lon.pop('ax2')
lon, lat = np.meshgrid(lon, lat)
aviso2['lat'] = lat
aviso2['lon'] = lon
aviso2['mlon'], aviso2['mlat'] = m2(aviso2['lon'], aviso2['lat'])
ssh = sp.loadmat('aviso_files/assh2.mat'); aviso2['ssh'] = ssh.pop('assh2')
u   = sp.loadmat('aviso_files/au2.mat'); aviso2['u'] = u.pop('au2')
v   = sp.loadmat('aviso_files/av2.mat'); aviso2['v'] = v.pop('av2')
# masking 
aviso2['u']   = apply_mask(aviso2['u'])
aviso2['v']   = apply_mask(aviso2['v'])
aviso2['ssh'] = apply_mask(aviso2['ssh'])
# scaling
aviso2['u'] = aviso2['u']/100; aviso2['v'] = aviso2['v']/100;

# HYCOM AREA 1 ========================================================
# Hardcoded for pathnames and filenames due to the nature of the 
# MAT files provided
hycom1 = {'lon':[],'lat':[],'mlon':[],'mlat':[],'u':[],'v':[],'ssh':[]}
lat = sp.loadmat('hycom_files/y1.mat'); lat = lat.pop('y1')
lon = sp.loadmat('hycom_files/x1.mat'); lon = lon.pop('x1')
lon, lat = np.meshgrid(lon, lat)
hycom1['lat'] = lat
hycom1['lon'] = lon
hycom1['mlon'], hycom1['mlat'] = m1(hycom1['lon'], hycom1['lat'])
u    = sp.loadmat('hycom_files/u1.mat'); hycom1['u'] = u.pop('u1')
v    = sp.loadmat('hycom_files/v1.mat'); hycom1['v'] = v.pop('v1')
ssh  = sp.loadmat('hycom_files/ssh1.mat'); hycom1['ssh'] = ssh.pop('ssh1')
# masking 
hycom1['u']   = apply_mask(hycom1['u'])
hycom1['v']   = apply_mask(hycom1['v'])
hycom1['ssh'] = apply_mask(hycom1['ssh'])
# scaling
hycom1['u'] = hycom1['u']/100; hycom1['v'] = hycom1['v']/100
hycom1['ssh'] = hycom1['ssh']*100 + 4

# HYCOM AREA 2 ========================================================
# Hardcoded for pathnames and filenames due to the nature of the 
# MAT files provided
hycom2 = {'lon':[],'lat':[],'mlon':[],'mlat':[],'u':[],'v':[],'ssh':[]}
lat = sp.loadmat('hycom_files/y2.mat'); lat = lat.pop('y2')
lon = sp.loadmat('hycom_files/x2.mat'); lon = lon.pop('x2')
lon, lat = np.meshgrid(lon, lat)
hycom2['lat'] = lat
hycom2['lon'] = lon
hycom2['mlon'], hycom2['mlat'] = m2(hycom2['lon'], hycom2['lat'])
u    = sp.loadmat('hycom_files/u2.mat'); hycom2['u'] = u.pop('u2')
v    = sp.loadmat('hycom_files/v2.mat'); hycom2['v'] = v.pop('v2')
ssh  = sp.loadmat('hycom_files/ssh2.mat'); hycom2['ssh'] = ssh.pop('ssh2')
# masking 
hycom2['u']   = apply_mask(hycom2['u'])
hycom2['v']   = apply_mask(hycom2['v'])
hycom2['ssh'] = apply_mask(hycom2['ssh'])
# scaling
hycom2['u'] = hycom2['u']/100; hycom2['v'] = hycom2['v']/100
hycom2['ssh'] = hycom2['ssh']*100



######################################################################
#    PLOTTING STUFF    ###############################################
######################################################################
print "\n\nPlotting %s\n\n" % ("."*50)

# FIGURE 1 ===========================================================
fig1 = plt.figure(facecolor='w',figsize=(8,12))
subplotfig(311, 1, 2.2, "aviso1", 1, "(a)", "AVISO")
subplotfig(312, 4, 2.2, "hycom1", 1, "(b)", "HYCOM")
subplotfig(313, 1, 2.2, "oe2", 1, "(c)", "OEII")
plt.savefig('leste2-vs-hycom-vs-aviso_2D.eps')
plt.savefig('leste2-vs-hycom-vs-aviso_2D.pdf')

# FIGURE 2 ===========================================================
fig2 = plt.figure(facecolor='w',figsize=(6,12))
subplotfig(311, 1, 3, "aviso2", 2, "(a)", "AVISO")
subplotfig(312, 5, 3, "hycom2", 2, "(b)", "HYCOM")
subplotfig(313, 2, 3, "proab", 2, "(c)", "PROAB")
plt.savefig('proab-vs-hycom-vs-aviso_2D.eps')
plt.savefig('proab-vs-hycom-vs-aviso_2D.pdf')

plt.show()

#######################################################################
#######################################################################

