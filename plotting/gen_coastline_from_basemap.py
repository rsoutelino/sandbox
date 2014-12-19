from mpl_toolkits.basemap import Basemap
import  matplotlib.pylab as plt



bbox = {'lat_min':-70,
                'lat_max':30,
                'lon_min':-107,
                'lon_max':19,}

m = Basemap(projection='cyl', llcrnrlat=bbox['lat_min'], urcrnrlat=bbox['lat_max'], 
            llcrnrlon=bbox['lon_min'], urcrnrlon=bbox['lon_max'], resolution='f')


coastline = m.drawcoastlines()

filename = open('coastline.txt','w')
filename.write('lon;lat;\n')
for p in range(len(coastline.axes.collections[0].get_paths())):
    path = coastline.axes.collections[0].get_paths()[p]
    vertices = path.vertices
    for i in vertices:
        filename.write('\n')
        for j in i:
            filename.write(str(j) + ';')
    filename.write('\n####################################################')
filename.close()




