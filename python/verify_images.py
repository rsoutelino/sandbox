#!/usr/bin/env python
######################################################
## Plots Blend SST for REMO Forecast
## Feb/2012, EDDIES Group - IEAPM
## rsoutelino@gmail.com
######################################################
import sys
import os
import re
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.io as sp
from   mpl_toolkits.basemap import Basemap

### DEFINING FUNCTIONS ###################################

def load_ascii(filename):
    SST = { 'lon':[],'lat':[],'mlon':[], 'mlat':[],'sst':[] }
    f = open(filename)
    lines = f.readlines()
    lon = lines[14]
    lat = lines[17]
    lon = lon.split(',')
    lat = lat.split(',')
    SST['lat'] = np.array(lat, dtype=float)
    SST['lon'] = np.array(lon, dtype=float)
    sst = lines[20:4021] # hardcoded for this size of sst array
    SST['sst'] = np.array( sst[0].split(',')[1:] , dtype=float )
    #print '\n\nReading SST array from %s %s\n\n' % (filename, '.'*20)
    for k in range(1, len(sst)):
        aux = np.array( sst[k].split(',')[1:] , dtype=float )
        SST['sst'] = np.vstack( (SST['sst'], aux) )
        
    SST['lon'], SST['lat'] = np.meshgrid(SST['lon'], SST['lat'])
    SST['sst'] =  SST['sst'] / 100 # scale factor
    SST['sst'] = pl.ma.masked_where(SST['sst'] < 0, SST['sst'])
    
    del f, lines, sst, lon, lat, aux
    
    return SST

def plotfig(size, filename, name, m1, m2, m3): 
    """
    size     : figure size [tuple]
    name     : dataset name [string]
    filename : filename of the dataset
    """
    etopo  = {'lon':[],'lat':[],'mlon':[],'mlat':[],'h':[], 'mh':[]}
    dataset = sp.loadmat( 'ETOPO_REMO_BR' )
    etopo['lon'] = dataset.pop('x')
    etopo['lat'] = dataset.pop('y')
    etopo['h']   = dataset.pop('z')
    etopo['lon'], etopo['lat'] = np.meshgrid(etopo['lon'], etopo['lat'])
    
    
    ### OPENING A FIGURE INSTANCE ################################################

    fig = plt.figure(facecolor='w', figsize=size)

    ### DEFINING THE AXES POSITIONS ##############################################

    ax1 = fig.add_axes([0.05, 0.05, 0.45, 0.9])
    ax2 = fig.add_axes([0.53, 0.515, 0.41, 0.42])
    ax3 = fig.add_axes([0.53, 0.07, 0.41, 0.42])
    ax4 = fig.add_axes([0.05, 0.046, 0.89, 0.015])
    
    ### PLOTTING SST AND ISOBATHS ################################################
    
    SST['mlon'], SST['mlat'] = m1(SST['lon'], SST['lat'])
    etopo['mlon'], etopo['mlat'] = m1(etopo['lon'], etopo['lat'])
    pc = m1.pcolormesh(SST['mlon'], SST['mlat'], SST['sst'], vmin=16, vmax=28, ax=ax1)
    m1.contour(etopo['mlon'], etopo['mlat'], etopo['h'], (-1000,-200),
               colors='w', linestyle='--', ax=ax1)
            
    SST['mlon'], SST['mlat'] = m2(SST['lon'], SST['lat'])
    etopo['mlon'], etopo['mlat'] = m2(etopo['lon'], etopo['lat'])
    pc = m2.pcolormesh(SST['mlon'], SST['mlat'], SST['sst'], vmin=16, vmax=28, ax=ax2)
    m2.contour(etopo['mlon'], etopo['mlat'], etopo['h'], (-1000,-200),
               colors='w', linestyle='--', ax=ax2)
            
    SST['mlon'], SST['mlat'] = m3(SST['lon'], SST['lat'])
    etopo['mlon'], etopo['mlat'] = m3(etopo['lon'], etopo['lat'])
    pc = m3.pcolormesh(SST['mlon'], SST['mlat'], SST['sst'], vmin=16, vmax=28, ax=ax3)
    m3.contour(etopo['mlon'], etopo['mlat'], etopo['h'], (-1000,-200),
               colors='w', linestyle='--', ax=ax3)
            
    ### CONTINENTS, COASTLINES, PARALLELS, MERIDIANS #############################
    
    m1.fillcontinents(color='k', ax=ax1)
    m2.fillcontinents(color='k', ax=ax2)
    m3.fillcontinents(color='k', ax=ax3)
    
    m1.drawparallels(np.arange(-28,-8,4),labels=[1, 0, 0, 0],
            dashes=[1,3], ax=ax1)
    m2.drawparallels(np.arange(-26,-19,1),labels=[0, 1, 0, 0],
            dashes=[1,3], ax=ax2)
    m3.drawparallels(np.arange(-25,-21,0.25),labels=[0, 1, 0, 0],
            dashes=[1,3], ax=ax3)
        
    m1.drawmeridians(np.arange(-48,-30,3),labels=[0, 0, 1, 0],
            dashes=[1,3], ax=ax1)
    m2.drawmeridians(np.arange(-48,-30,2),labels=[0, 0, 1, 0],
            dashes=[1,3], ax=ax2)
    m3.drawmeridians(np.arange(-48,-30,0.5),labels=[0, 0, 1, 0],
            dashes=[1,3], ax=ax3)

    ### BOUNDING BOXES ###########################################################

    lx, ly = m1( (llc2x, llc2x, urc2x, urc2x, llc2x),
                 (llc2y, urc2y, urc2y, llc2y, llc2y) )
    m1.plot(lx, ly, color='0.4', linestyle='--', ax=ax1)

    lx, ly = m2( (llc3x, llc3x, urc3x, urc3x, llc3x),
                 (llc3y, urc3y, urc3y, llc3y, llc3y) )
    m2.plot(lx, ly, color='0.4', linestyle='--', ax=ax2)
    
    ### TITLE ####################################################################
    
    tx, ty = m1(-48.4,-11)
    ax1.text(tx, ty,'GHRSST 1km Blended SST - G1SST', color='w',
            style='italic')
    tx, ty = m1(-45, -11.7)
    txt = filename[10:12]+'/'+filename[8:10]+'/'+filename[4:8]
    ax1.text(tx, ty, txt, color='w')

    ### COLORBAR ##############################################################      
  
    cbar = plt.colorbar(pc,cax=ax4,orientation='horizontal')
    cbar.set_label('$^\circ$ C')

    return filename






### SCRIPT STARTS HERE ########################################################

### DEFINING MAP PROJECTIONS #################################################

print "\n\n DEFINING MAP PROJECTIONS %s %s" % (file, "."*50)

llc1y = -30; urc1y = -10
llc2y = -26; urc2y = -20
llc3y = -23.5; urc3y = -22.2
llc1x = -50; urc1x = -35
llc2x = -48; urc2x = -39
llc3x = -43.2; urc3x = -41.3

m1 = Basemap(projection='merc',llcrnrlat=llc1y,urcrnrlat=urc1y,
             llcrnrlon=llc1x,urcrnrlon=urc1x,lat_ts=0,resolution='i')
                
m2 = Basemap(projection='merc',llcrnrlat=llc2y,urcrnrlat=urc2y,
             llcrnrlon=llc2x,urcrnrlon=urc2x,lat_ts=0,resolution='h')
                
m3 = Basemap(projection='merc',llcrnrlat=llc3y,urcrnrlat=urc3y,
             llcrnrlon=llc3x,urcrnrlon=urc3x,lat_ts=0,resolution='f')

base_dir = os.path.abspath(os.path.dirname(__file__))

file = open(os.path.join(base_dir, 'sst_20090401.png',), 'r')

for file in os.listdir(base_dir):
    if re.search(r'\.(ascii)$',file):
        prefix = re.split(r'\..',file)[0]
        png = prefix + '.png'
        if png not in os.listdir(base_dir):
            print "\n\n LOADING %s %s" % (file, "."*50)
            SST = load_ascii(file)
            print "Plotting maps"
            plotfig( (13, 10), file, SST, m1, m2, m3 )
            print "Saving PNG figure"
            plt.savefig("%s.png" %(file[:12]) )
            plt.close('all')



