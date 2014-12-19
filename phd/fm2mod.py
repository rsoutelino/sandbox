#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script for writing modfile to HOPS OA
#
# Rafael Soutelino - rsoutelino@gmail.com
#
#
# Last modification: Apr, 2011
#####################################################################

print ' \n' + '==> ' + '  IMPORTING MODULES ...\n' + ' ' 

# IMPORTING MODULES #################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
# import datetime as dt
import netCDF4 as nc
import scipy.io as sp
from datetime import datetime

def near(x,x0):
        """
        Find the index where x has the closer value to x0
        """

        dx = x - x0
        dx = np.abs(dx)
        fn = np.where( dx == dx.min() )
        fn = fn[0][0]

        return fn

#####################################################################

# filename
f = open('bc-nbuc_fm.mod', 'w')

# loading data
dataset = sp.loadmat('BC-NBUC_FM.mat')
lon   = dataset['lon'][0,:]
lat   = dataset['lat'][:,0]
depth = dataset['z'][:]; depth = np.abs(np.squeeze(depth))

# getting a subdomain
fj1 = near(lon, -42); fj2 = near(lon, -30);
fi1 = near(lat, -23); fi2 = near(lat, -5);

# re-loading data
d = 1 # sub-sampling factor
lon   = dataset['lon'][fi1:fi2:d,fj1:fj2:d]
lat   = dataset['lat'][fi1:fi2:d,fj1:fj2:d]
U     = dataset['u'][:,fi1:fi2:d,fj1:fj2:d]
V     = dataset['v'][:,fi1:fi2:d,fj1:fj2:d]
temp  = dataset['temp'][:,fi1:fi2:d,fj1:fj2:d] 
salt  = dataset['salt'][:,fi1:fi2:d,fj1:fj2:d]

li, lj = lon.shape

# mod file header
line  = ' title = BC-NBUC Feature Model';                 f.write(line + '\n')
line  = ' stations = ' + str(li * lj);                    f.write(line + '\n')
line  = ' str_time = 14833.0000, Jan 01 2009 00:00:00';   f.write(line + '\n')
line  = ' end_time = 14833.0000, Jan 01 2009 00:00:00';   f.write(line + '\n')
line  = ' Jday_offset = 2440000';                         f.write(line + '\n')
line  = ' lng_min = ' + str(lon.min());                   f.write(line + '\n')                           
line  = ' lng_max = ' + str(lon.max());                   f.write(line + '\n')           
line  = ' lat_min = ' + str(lat.min());                   f.write(line + '\n')                          
line  = ' lat_max = ' + str(lat.max());                   f.write(line + '\n')
line  = ' format = ascii, record interleaving';           f.write(line + '\n')      
line  = ' type = CTD';                                    f.write(line + '\n')       
line  = ' fields = depth, temperature, salinity';         f.write(line + '\n')     
line  = ' units = meter, Celcius, PSU';                   f.write(line + '\n')   
line  = ' creation_date = Tue Apr 19 10:11:20 2011';      f.write(line + '\n')
line  = 'END';                                            f.write(line + '\n')


# writing data arrays into mod file
#long, latg = np.meshgrid(lon, lat)
long = lon; latg = lat
c = 0

for i in range(0, li-1):
    for j in range(0, lj-1):
        t = temp[:,i,j]
        s = salt[:,i,j]
        F = np.where(s < 40)
        if np.size(F) == 0:
            pass
        else: 
            c = c + 1
            z = depth[F[0][0]:F[0][-1]+1]; z = np.round(z*10)
            t = t[F[0][0]:F[0][-1]+1];  t = np.round(t*1000)
            s = s[F[0][0]:F[0][-1]+1];  s = np.round(s*1000)
            
            f.write(' 3    ' + str(z.size) +'     '+ str(c) +'  '+ str(long[i,j]) +' '+ str(latg[i,j])
                    +'   '+ str(z.max()/10) + '  14833.0000 1.00e-01 1.00e-03 1.00e-03\n')
            f.write("  0 'CTD: z t s' \n")
            
            S = t.size
            
            if S <= 10:
                z = np.array(z,dtype=int); z1 = str(z)
                t = np.array(t,dtype=int); t1 = str(t)
                s = np.array(s,dtype=int); s1 = str(s)          
                f.write(z1[1:-1] + '\n')
                f.write(t1[1:-1] + '\n')
                f.write(s1[1:-1] + '\n') 
                
            else:
            
                lines = S / 10
                lastline = np.remainder(S,10)
                if lastline != 0:
                    z1 = z[:lastline*-1]; z2 = z[lastline*-1:]; z2 = np.array(z2,dtype=int)
                    t1 = t[:lastline*-1]; t2 = t[lastline*-1:]; t2 = np.array(t2,dtype=int)
                    s1 = s[:lastline*-1]; s2 = s[lastline*-1:]; s2 = np.array(s2,dtype=int)
                    z1.shape = (lines, 10); z1 = np.array(z1,dtype=int)
                    t1.shape = (lines, 10); t1 = np.array(t1,dtype=int)
                    s1.shape = (lines, 10); s1 = np.array(s1,dtype=int)
                    
                    for l in range(0, lines):
                        z11 = str(z1[l,:]); f.write(z11[1:-1] + '\n')
                    
                    z21 = str(z2); f.write(z21[1:-1] + '\n')
                       
                    for l in range(0, lines):
                        t11 = str(t1[l,:]); f.write(t11[1:-1] + '\n')
                        
                    t21 = str(t2); f.write(t21[1:-1] + '\n')
                        
                    for l in range(0, lines):
                        s11 = str(s1[l,:]); f.write(s11[1:-1] + '\n')
                        
                    s21 = str(s2); f.write(s21[1:-1] + '\n')
                                                           
                else:
                    z.shape = (lines, 10); z = np.array(z,dtype=int)
                    t.shape = (lines, 10); t = np.array(t,dtype=int)
                    s.shape = (lines, 10); s = np.array(s,dtype=int)
                
                    for l in range(0, lines):
                        z1 = str(z[l,:]); f.write(z1[1:-1] + '\n')
                    
                    for l in range(0, lines):
                        t1 = str(t[l,:]); f.write(t1[1:-1] + '\n')
                    
                    for l in range(0, lines):
                        s1 = str(s[l,:]); f.write(s1[1:-1] + '\n') 
                    
                      
                        
print "PLEASE REMEMBER TO CHANGE THE NUMBER OF STATIONS IN THE HEADER TO: " + str(c)

print "I STILL NEED TO IMPLEMENT THIS BETTER BY LEARN HOW TO OVERWRITE THINGS TO"
print "PRE_EXISTING FILES"













