######################################################
## plotando os dados do ADCP processados pelo CODAS ##
##         Rafael Soutelino - set/2007              ##
######################################################
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import scipy.io as sp
from mpl_toolkits.basemap import Basemap

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
    
#################################################
# preparing to plot figures

m = Basemap(projection='merc',llcrnrlat=-22.5,urcrnrlat=-12.0,llcrnrlon=-41.0,urcrnrlon=-33.7,lat_ts=0,resolution='l')
s = 2; # escala para vetores
ct = 200; # ct2 = 20;
pp = 1000;
p = 50
etopo = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/etopo2_leste.mat')
xb = etopo.pop('xb')
yb = etopo.pop('yb')
zb = etopo.pop('zb')
xb,yb = m(xb,yb)
# masking topography to plot continental shelf mask
zb = pl.ma.masked_where(zb <= -pp, zb)


#################################################
fig1 = plt.figure(1,figsize=(15,8),facecolor='w')

###################################################################################

oe1 = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/psigeo_OEI.mat')
psigo = oe1.pop('psigo')
xgi = oe1.pop('xgi'); ygi = oe1.pop('ygi');
xgi = xgi[:-10, :]; ygi = ygi[:-10, :];
x, y = m(xgi, ygi)

dataset = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/occam_2003.mat')
lon = dataset.pop('lon'); lat = dataset.pop('lat');
u   = dataset.pop('u')  ; v   = dataset.pop('v')  ;

u = griddata(lon.ravel(), lat.ravel(), u.ravel(), x, y)
v = griddata(lon.ravel(), lat.ravel(), v.ravel(), x, y)

p1 = plt.subplot(131)

m.quiver(x, y, u*s, v*s,scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[1, 0, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p1.text(ax, ay, "2003 Mean", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p1.text(ax, ay, "OCCAM", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p1.text(ax, ay, "a", ha="center", va="center", size=12,bbox=bbox_props)

######################################################################

dataset = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/occam_feb2003.mat')
lon = dataset.pop('lon'); lat = dataset.pop('lat');
u   = dataset.pop('u')  ; v   = dataset.pop('v')  ;

u = griddata(lon.ravel(), lat.ravel(), u.ravel(), x, y)
v = griddata(lon.ravel(), lat.ravel(), v.ravel(), x, y)

p2 = plt.subplot(132)

m.quiver(x, y, u*s, v*s,scale=10);
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
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p2.text(ax, ay, "Feb 2003", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p2.text(ax, ay, "OCCAM", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p2.text(ax, ay, "b", ha="center", va="center", size=12,bbox=bbox_props)

#############################################################################################

dataset = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/occam_sep2003.mat')
lon = dataset.pop('lon'); lat = dataset.pop('lat');
u   = dataset.pop('u')  ; v   = dataset.pop('v')  ;

u = griddata(lon.ravel(), lat.ravel(), u.ravel(), x, y)
v = griddata(lon.ravel(), lat.ravel(), v.ravel(), x, y)

p3 = plt.subplot(133)

m.quiver(x, y, u*s, v*s,scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
m.drawparallels(np.arange(-22,-12,2),labels=[0, 1, 0, 0],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
text = np.str(p) + ' m'
ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p3.text(ax, ay, "Sep 2003", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p3.text(ax, ay, "OCCAM", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p3.text(ax, ay, "c", ha="center", va="center", size=12,bbox=bbox_props)

#############################################################################################
fig1.tight_layout(pad=3)
plt.show()
plt.savefig('figures/occam.pdf')
# plt.savefig('figures/occam.eps')
  


