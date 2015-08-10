#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import time

dt = 300
plt.figure()

os.system('sh roms_ek-1.sh')
data = np.loadtxt('phd16_ek-1.out')
t = data[:,0]
t = (t*dt)  / (86400*360 )
ek = data[:,5]
plt.plot(t, ek, 'b')
plt.grid('on')
plt.title('Total KE temporal evolution')
plt.xlabel('Years of simulation')
plt.ylabel('J m$^{-2}$')

os.system('sh roms_ek-2.sh')
data2 = np.loadtxt('phd16_ek-2.out')
t2  = data2[:,0]
ek2 = data2[:,5]
f = np.where(t2 >= 400000)
t2[f], ek2[f] = np.nan, np.nan
t2 = ((t2*dt)  / (86400*360 )) + t[-1]
plt.plot(t2, ek2, 'r')

os.system('sh roms_ek-3.sh')
data3 = np.loadtxt('phd16_ek-3.out')
t3  = data3[:,0]
ek3 = data3[:,5]
f = np.where(t3 >= 400000)
t3[f], ek3[f] = np.nan, np.nan
t3 = ((t3*dt)  / (86400*360 )) + t2[-1]
plt.plot(t3, ek3, 'g')

plt.show()


	

