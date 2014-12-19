# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
import Image
from cookb_signalsmooth import smooth


def draw_logo(fig, logo_position, filename):
        """
        fig : matplotlib figure object
        axes: list with axes position and size
        """              
        ax2 = fig.add_axes(logo_position)
        logofile = Image.open(filename)
        logo = np.asarray(logofile)
        im = plt.imshow(logo)
        ax2.set_axis_off()


################################################################################

vel = sp.loadmat('veloc_barocl_xbt.mat')
var = sp.loadmat('xbt_variables.mat')

xif, zbf, lats, lons = (var['xif'], var['zbf'], var['lats'], var['lons'])

z = vel['z']
lat2 = vel['lat2'].ravel()
lon2 = vel['lon2'].ravel()
vg = vel['velg']

lat = []
lon = []

for k in range(lat2.shape[-1] - 1):
    lat.append(lat2[k:k+2].mean())
    lon.append(lon2[k:k+2].mean())


lat, z = np.meshgrid(lat, z)
lon = np.array(lon)


fig4 = plt.figure(facecolor='w', figsize=(14,8))
plt.contourf(lat, z, vg, np.arange(-0.8, 0.8+0.06, 0.06), cmap=plt.cm.RdBu)
cbar = plt.colorbar()
plt.fill(xif, zbf, '0.7', zorder=3)
plt.plot(lats, 0*lats, 'wv', markersize=7)
plt.axis([lats[0]+0.01, lats[-1]-0.01, 1200, -20])
cbar.set_label('m s$^{-1}$')
plt.xlabel('Latitude [$^o$S]')
plt.ylabel('Profundidade [m]')
plt.title(u'Velocidade Geostrófica: UH-14 22/09/2012')
draw_logo(fig4, [0.12, 0.14, 0.15, 0.15], 'logo-left.png')
draw_logo(fig4, [0.23, 0.14, 0.15, 0.15], 'photos/IMG_5274.JPG')
plt.savefig('eddy_hunter_gvel.png', dpi=300)
plt.show()




area = 1200 * np.abs( ( lats[0] - lats[-1] ) * 111000 )[0]
vm   = vg.mean()

tv_total = (area * vm) / 1e6

vg_cb = np.ma.masked_where(vg > 0, vg)
vm_cb = vg_cb.mean()

tv_cb = ((area/2) * vm_cb) / 1e6


lat = lat[0,:]
v = vg[0,:]


v, lon, lat, lons, lats = v[::-1], lon[::-1], lat[::-1], lons[::-1], lats[::-1]


x = np.linspace(lons[0], lon[-1], 15)
y = np.linspace(lats[0], lat[-1], 15)



vi = np.interp(y, lat, v)
vi = smooth(vi, window_len=6, window='hanning')

u = vi.copy()
v = u*0

mdict = {'u':u, 'v':v, 'x':x, 'y':y}
sp.savemat('gvel_map.mat', mdict)