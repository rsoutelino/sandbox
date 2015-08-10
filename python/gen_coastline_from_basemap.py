from mpl_toolkits.basemap import Basemap
import  matplotlib.pylab as plt
import scipy.io as sp
import numpy as np
import xray


bbox = {'lat_min':-48,
                'lat_max':-33,
                'lon_min':165,
                'lon_max':180,}

m = Basemap(projection='cyl', llcrnrlat=bbox['lat_min'], urcrnrlat=bbox['lat_max'], 
            llcrnrlon=bbox['lon_min'], urcrnrlon=bbox['lon_max'], resolution='h')


coastline = m.drawcoastlines()

lon, lat = [], []
filename = open('coastline.txt','w')
filename.write('lon;lat;\n')
for p in range(len(coastline.axes.collections[0].get_paths())):
    path = coastline.axes.collections[0].get_paths()[p]
    vertices = path.vertices
    for i in vertices:
        filename.write('\n')
        for j in i:
            filename.write(str(j) + ';')
            lon.append(i[0])
            lat.append(i[1])
    filename.write('\n####################################################')
filename.close()

vdict = dict(lon=lon, lat=lat)
coastline = xray.Dataset(vdict)
coastline.to_netcdf('coastline.nc')


