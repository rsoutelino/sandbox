#!/usr/local/epd-7.2-1-rh5-x86_64/bin/python
# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
from romslab import subset, LoadEtopo5
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import seawater.csiro as sw
import Image
from cookb_signalsmooth import smooth
################################################################################

def mask_nans(array, missing_value=-99999):
    f = np.where(np.isnan(array)==1)
    array[f] = missing_value
    array = np.ma.masked_where(array==missing_value, array)
    return array

def refine(lon, lat, sst, lims, res=0.02):
	xg = np.arange(lims[0]-1, lims[1]+1, res)
	yg = np.arange(lims[2]-1, lims[3]+1, res)
	xg, yg = np.meshgrid(xg, yg)
	sst = griddata(lon.ravel(), lat.ravel(), sst.ravel(), xg, yg)
	return xg, yg, sst

def draw_logo(fig, logo_position):
        """
        fig : matplotlib figure object
        axes: list with axes position and size
        """              
        ax2 = fig.add_axes(logo_position)
        logofile = Image.open('logo-left.png')
        logo = np.asarray(logofile)
        im = plt.imshow(logo)
        ax2.set_axis_off()

################################################################################
xmin, xmax, ymin, ymax = -44, -40.5, -25, -22
lims = (xmin, xmax, ymin, ymax)
m = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,
             llcrnrlon=xmin,urcrnrlon=xmax,lat_ts=0,resolution='h')

topo = sp.loadmat('/home/rsoutelino/misc/etopo/etopo1_atlsul.mat')
xb, yb, zb = topo['lon'], topo['lat'], topo['topo']
xb, yb = np.meshgrid(xb, yb)
xb, yb, zb = subset(xb, yb, zb, xmin, xmax, ymin, ymax)

filename = 'mosaico_20120919.mat'
data = sp.loadmat(filename)
lon, lat, sst = data['lon'], data['lat'], data['sst']
sst = mask_nans(sst, missing_value=18)
f = np.where(sst < 0)
sst[f] = 18
lon, lat, sst = refine(lon, lat, sst, lims, res=0.02)

base = -42.097463607788086, -22.817603940087874 
outA = -42.296874999999993, -24.568749999999998
outB = -42.66718749999999, -24.245312499999997
inA  = -42.146874999999994, -23.237499999999997
inB  = -42.254687499999996, -23.209374999999994
eddy = np.loadtxt('eddy.txt')
eddy[:,0] = smooth(eddy[:,0], window_len=11, window='hanning')
eddy[:,1] = smooth(eddy[:,1], window_len=11, window='hanning')

trajectoryA = [np.linspace(base[0], outA[0], 10), 
              np.linspace(base[1], outA[1], 10)]
trajectoryB = [np.linspace(base[0], outB[0], 10), 
              np.linspace(base[1], outB[1], 10)]

stationsA = [np.linspace(inA[0], outA[0], 6), np.linspace(inA[1], outA[1], 6)]
stationsB = [np.linspace(inB[0], outB[0], 6), np.linspace(inB[1], outB[1], 6)]

distA, distB = [], []

for k in range(stationsA[0].size):
	distA.append( int(sw.dist([base[0], stationsA[0][k]], [base[1], stationsA[1][k]], 'nm')[0]) )
	distB.append( int(sw.dist([base[0], stationsB[0][k]], [base[1], stationsB[1][k]], 'nm')[0]) )

fig = plt.figure(facecolor='w', figsize=(15,10))
pc = m.pcolormesh(lon, lat, sst, vmin=18, vmax=24, cmap=plt.cm.Spectral_r)
# pc = m.contourf(lon, lat, sst, np.arange(18, 25, 0.1), cmap=plt.cm.Spectral_r)
m.contour(xb, yb, zb, (-1000, -200), alpha=0.6)
m.fillcontinents()
m.drawcoastlines()
m.drawparallels(np.arange(-25,-21,0.5),labels=[1, 0, 0, 0],
            dashes=[1,1000])
m.drawmeridians(np.arange(-48,-30,1),labels=[1, 0, 1, 0],
            dashes=[1,1000])

trajectory, stations, dist = trajectoryA, stationsA, distA

# m.plot(trajectory[0], trajectory[1], 'k')
m.plot(stations[0], stations[1], 'wv', markersize=8)
m.plot(base[0], base[1], '*w', markersize=10)
# bbox_props = dict(boxstyle="round", fc="k", ec="0.5", alpha=0.2)
# for k in range(stations[0].size):
# 	plt.text(stations[0][k]+0.05, stations[1][k]-0.02, "%s nm" %dist[k],
# 	         va="center", size=10, bbox=bbox_props)
m.plot(eddy[:,0], eddy[:,1], 'k', alpha=0.2, linewidth=10)

gvel = sp.loadmat('gvel_map.mat')
m.quiver(gvel['x'], gvel['y'], gvel['u'], gvel['v'], scale=3, alpha=0.5)
         # width=0.01, headwidth=2, alpha=0.5)


ax1 = fig.add_axes([0.34, 0.85, 0.25, 0.015])
cbar = plt.colorbar(pc, cax=ax1, orientation='horizontal', extend='both')
plt.title(filename[:-4] + u" [$^o$C]", fontsize=11)
draw_logo(fig, [0.20, 0.75, 0.12, 0.12])
# m.plot(lons, lats, 'w*')
real_stations = np.loadtxt('positions_real.dat')
# m.plot(real_stations[:,1], real_stations[:,2], '*b', markersize=8)


plt.show()
plt.savefig("gvel_map.png", dpi=300)