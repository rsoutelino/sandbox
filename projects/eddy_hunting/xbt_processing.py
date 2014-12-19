# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
from romslab import near
from matplotlib.mlab import griddata
from matplotlib import delaunay
import Image
import glob
from cookb_signalsmooth import smooth
################################################################################

class ReadXBT(object):
	"""docstring for ReadXBT"""
	def __init__(self, filename):
		self.filename = filename
		f = open(filename)
		data = f.readlines()
		data = data[36:]
		self.id = int(self.filename[11])
		self.temp, self.pres = [], []
		for line in data:
			if 'station1' in self.filename:
				self.pres.append( float(line.split('\t')[2]) )
				self.temp.append( float(line.split('\t')[3]) )
			else:
				self.pres.append( float(line.split('\t')[0]) )
				self.temp.append( float(line.split('\t')[1]) )
		self.pres = np.array(self.pres, dtype=np.float64)
		self.temp = np.array(self.temp, dtype=np.float64)
        

def remove_nans(array):
	f = np.where(np.isnan(array) == 0)
	array = array[f[0]]
	return array

def filter(temp, pres):
	temp, pres = remove_nans(temp), remove_nans(pres)
	f = np.where(temp > 2)
	temp, pres = temp[f], pres[f]
	k = 0
	a = 5
	b = 0.22
	for j in range(temp.size):
		try:
			if temp[k] > temp[k-a:k+a].mean() + b or temp[k] < temp[k-a:k+a].mean() - b:
				temp[k], pres[k] = np.nan, np.nan
				temp, pres = remove_nans(temp), remove_nans(pres)
		except IndexError:
			pass
		k += 1
	return temp, pres


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

filelist = glob.glob('xbt/*.edf')
filelist.sort()
filelist.remove('xbt/station2_T10.edf')
filelist.remove('xbt/station1_T5.edf')
xbt = []
for station in filelist:
	xbt.append( ReadXBT(station) )
pos = np.loadtxt('positions_real.dat')
stations, lons, lats, pmax = pos[:-1,0], pos[:-1,1], pos[:-1,2], pos[:-1,3]


fig1 = plt.figure(facecolor='w', figsize=(10,10))
for k in range(len(xbt)):
	plt.plot(xbt[k].temp, -xbt[k].pres, label=xbt[k].filename[4:16])
plt.legend(loc=4)
plt.grid()
plt.title(u'Estações de XBT Brutas')
plt.savefig('brutos.png', dpi=300)
plt.show()

fig2 = plt.figure(facecolor='w', figsize=(10,10))
for k in range(len(xbt)):
	if 'station1_T5' in xbt[k].filename:
		# xbt[k].temp, xbt[k].pres = filter(xbt[k].temp, xbt[k].pres)
		# xbt[k].temp, xbt[k].pres = filter(xbt[k].temp, xbt[k].pres)
		xbt[k].temp = smooth(xbt[k].temp, window_len=201, window='hanning')

	elif 'station2_T5' in xbt[k].filename:
		xbt[k].temp, xbt[k].pres = filter(xbt[k].temp, xbt[k].pres)
		xbt[k].temp, xbt[k].pres = filter(xbt[k].temp, xbt[k].pres)
		xbt[k].temp = smooth(xbt[k].temp, window_len=51, window='hanning')
	else:
		xbt[k].temp, xbt[k].pres = filter(xbt[k].temp, xbt[k].pres)
		xbt[k].temp = smooth(xbt[k].temp, window_len=31, window='hanning')

	plt.plot(xbt[k].temp, -xbt[k].pres, label=xbt[k].filename[4:16])
plt.legend(loc=4)
plt.grid()
plt.title(u'Estações de XBT Filtradas')
plt.savefig('filtrados.png', dpi=300)
plt.show()


z = np.arange(0,1800)
xx, zz = np.meshgrid(lats, z)
k = 2

flag = []

for st in xbt:
	st.temp = np.interp(z, st.pres, st.temp)
	f1 = np.where(stations == st.id)[0]
	flag.append(z < pmax[f1])
	st.temp.shape = (st.temp.size, 1)
	k += 1 

flag = np.array(flag)
flag = flag.transpose()
flag = np.fliplr(flag)
Tsec = np.concatenate( (xbt[3].temp, xbt[2].temp, xbt[1].temp, xbt[0].temp), axis=1 )

Tsec = np.ma.masked_where(flag == False, Tsec)
xx   = np.ma.masked_where(flag == False, xx)
zz   = np.ma.masked_where(flag == False, zz)
zi, xi = z[::5], np.linspace(lats[0], lats[-1], 200)
xg, zg = np.meshgrid(xi, zi)

tri = delaunay.Triangulation(xx.compressed(), zz.compressed())
interpolate = tri.nn_extrapolator(Tsec.compressed())
Tg = interpolate(xg, zg)

zb = np.interp(xi[::-1], lats[::-1], pmax[::-1])
zb = zb[::-1]
zb = smooth(zb, window_len=51, window='hanning')

xif = np.hstack( ( xi[0], xi, xi[-1], xi[0] ) )
zbf = np.hstack( ( zb[0], zb, zb.max()+100, zb.max()+100 ) )

fig3 = plt.figure(facecolor='w', figsize=(14,8))
plt.contourf(xg, zg, Tg, 40, cmap=plt.cm.Spectral_r)
cbar = plt.colorbar()
plt.contour(xg, zg, Tg, 20, colors='k', zorder=1)
plt.fill(xif, zbf, '0.7', zorder=3)
plt.plot(lats, 0*lats, 'wv', markersize=7)
plt.axis([lats[0]+0.01, lats[-1]-0.01, 1200, -20])
# plt.axis([lats[0]+0.01, -24.445625, 1200, -20])
# plt.axis([lats[0]+0.01, -24.56875, 1200, -20])
cbar.set_label('$^o$C')
plt.xlabel('Latitude [$^o$S]')
plt.ylabel('Profundidade [m]')
plt.title('XBT: UH-14 22/09/2012')
draw_logo(fig3, [0.12, 0.14, 0.15, 0.15], 'logo-left.png')
draw_logo(fig3, [0.23, 0.14, 0.15, 0.15], 'photos/IMG_5274.JPG')
plt.savefig('eddy_hunter_temp.png', dpi=300)
plt.show()

mdict = {'xif':xif, 'zbf':zbf, 'lats':lats, 'lons':lons}
sp.savemat('xbt_variables.mat', mdict)

