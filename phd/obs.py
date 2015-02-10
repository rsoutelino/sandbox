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
# lon,lat = m(lon,lat)
# masking topography to plot continental shelf mask\
zb = pl.ma.masked_where(zb <= -pp, zb)


#################################################
fig1 = plt.figure(1,figsize=(15,8),facecolor='w')

###################################################################################

oe1 = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/psigeo_OEI.mat')
psigo = oe1.pop('psigo')
xgi = oe1.pop('xgi'); ygi = oe1.pop('ygi');
xgi,ygi = m(xgi,ygi)
uctd = oe1.pop('uctd'); vctd = oe1.pop('vctd');  

# masking 
f = pl.find(pl.isnan(psigo) == 1); 
a, b = psigo.shape; 
psigo = pl.ravel(psigo); uctd = pl.ravel(uctd); vctd = pl.ravel(vctd);

psigo[f] = 0; uctd[f] = 0; vctd[f] = 0;
psigo.shape = (a,b); uctd.shape = (a,b); vctd.shape = (a,b)
psigo = pl.ma.masked_where(psigo == 0, psigo)
uctd = pl.ma.masked_where(uctd == 0, uctd)
vctd = pl.ma.masked_where(vctd == 0, vctd)

g = 1

p1 = plt.subplot(131)
#con = m.contourf(xgi,ygi,psigo,scale,linewidth=0,alpha = 0.5); 
m.quiver(xgi[::g,::g],ygi[::g,::g],uctd[::g,::g]*s,vctd[::g,::g]*s,scale=10);
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
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p1.text(ax, ay, "Dec 2001", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p1.text(ax, ay, "T,S-derived", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p1.text(ax, ay, "a", ha="center", va="center", size=12,bbox=bbox_props)

######################################################################

oe2 = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/psiobs_OEII.mat')
psiob = oe2.pop('psiob')
xgi = oe2.pop('xgi'); ygi = oe2.pop('ygi');
xgi,ygi = m(xgi,ygi)
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

p2 = plt.subplot(132)
#m.contourf(xgi,ygi,psiob,scale,linewidth=0,alpha = 0.5);
m.quiver(xgi[::g,::g],ygi[::g,::g],uo[::g,::g]*s,vo[::g,::g]*s,scale=10);
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
  #[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
  #set(h,'facecolor',[.7 .7 .7]);
  #[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
  #set(h,'color',[.7 .7 .7]);
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
p2.text(ax, ay, "Feb 2005", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-37.2,-22);
p2.text(ax, ay, "ADCP-derived", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p2.text(ax, ay, "b", ha="center", va="center", size=12,bbox=bbox_props)

#############################################################################################

proab1 = sp.loadmat('/home/rsoutelino/prod/grl_msc/proc/psiobs_proab.mat')
psiob = proab1.pop('psiob')
xgi = proab1.pop('xgi'); ygi = proab1.pop('ygi');
xgi,ygi = m(xgi,ygi)
uo = proab1.pop('uo'); vo = proab1.pop('vo');  

# loading drifters information
d16 = pl.loadtxt('/home/rsoutelino/projetos/proabrolhos/derivadores/bd69016.txt')
d20 = pl.loadtxt('/home/rsoutelino/projetos/proabrolhos/derivadores/bd69020.txt')
lat16 = d16[:,1]; lon16 = d16[:,2]; yy16 = d16[:,4]; mm16 = d16[:,5]; dd16 = d16[:,7]
lat20 = d20[:,1]; lon20 = d20[:,2]; yy20 = d20[:,4]; mm20 = d20[:,5]; dd20 = d20[:,7]
lon16,lat16 = m(lon16,lat16)
lon20,lat20 = m(lon20,lat20)

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

p3 = plt.subplot(133)
#m.contourf(xgi,ygi,psiob,scale,linewidth=0,alpha = 0.5);
m.quiver(xgi[::g,::g],ygi[::g,::g],uo[::g,::g]*s,vo[::g,::g]*s,scale=10);
m.plot(lon16,lat16,color='k',linewidth=2)
m.plot(lon20,lat20,color='k',linewidth=2)
m.contourf(xb,yb,zb,colors=('0.9'),linewidth=0);
m.drawcoastlines(zorder=5)
m.drawcountries(zorder=4)
m.fillcontinents(color='0.8',lake_color='aqua',zorder=3)
  #[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
  #set(h,'facecolor',[.7 .7 .7]);
  #[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
  #set(h,'color',[.7 .7 .7]);
m.drawparallels(np.arange(-22,-12,2),labels=[0, 0, 0, 1],dashes=[1,1000],zorder=6)
m.drawmeridians(np.arange(-41,-34,2),labels=[0, 0, 1, 1],dashes=[1,1000],zorder=7)
text = "50 cm/s";
ax,ay = m(-40.5,-13.0); m.quiver(ax,ay,0.5*s,0,scale=10,color='k',zorder=10)
ax,ay = m(-40.5,-13.4); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
#text = np.str(p) + ' m'
#ax,ay = m(-34.8,-22.0); plt.text(ax,ay,text,color='k',fontsize=10,fontweight='bold')
ax,ay = m(-39.5,-18.5); plt.text(ax,ay,"ABROLHOS",color='k',fontsize=10)
ax,ay = m(-39.2,-18.8); plt.text(ax,ay,"BANK",color='k',fontsize=10)
ax,ay = m(-38.8,-16.2); plt.text(ax,ay,"RCB",color='k',fontsize=10)
ax,ay = m(-40,-14.5); plt.text(ax,ay,"Ilheus",color='k',fontsize=10)
ax,ay = m(-39.3,-12.5); plt.text(ax,ay,"Salvador",color='k',fontsize=10)
ax,ay = m(-40.75,-17.6); plt.text(ax,ay,"Caravelas",color='k',fontsize=10)
ax,ay = m(-38,-21); plt.text(ax,ay,"VITORIA-TRINDADE RIDGE",color='k',fontsize=10)
bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
ax,ay = m(-34.7,-12.5);
p3.text(ax, ay, "Sep 2007", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-36,-16.5);
p3.text(ax, ay, "ADCP-derived at 50 m", ha="center", va="center", size=10,bbox=bbox_props)
ax,ay = m(-36.5,-22);
p3.text(ax, ay, "Surface drifters", ha="center", va="center", size=10,bbox=bbox_props)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
ax,ay = m(-40.5,-22)
p3.text(ax, ay, "c", ha="center", va="center", size=12,bbox=bbox_props)

#############################################################################################
fig1.tight_layout(pad=3)
plt.show()
plt.savefig('figures/obs.pdf')
# plt.savefig('figures/obs.eps')

  


