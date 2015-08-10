#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import time

dt = 300


data = np.loadtxt('phd15_ek-1.out')
t = data[:,0]
t = (t*dt)  / (86400*360 )
ek = data[:,5]

data2 = np.loadtxt('phd15_ek-2.out')
t2  = data2[:,0]
ek2 = data2[:,5]
f = np.where(t2 >= 400000)
t2[f], ek2[f] = np.nan, np.nan
t2 = ((t2*dt)  / (86400*360 )) + t[-1]

data3 = np.loadtxt('phd15_ek-3.out')
t3  = data3[:,0]
ek3 = data3[:,5]
f = np.where(t3 >= 400000)
t3[f], ek3[f] = np.nan, np.nan
t3 = ((t3*dt)  / (86400*360 )) + t2[-1]

t_phd15 = np.concatenate( ( t, t2, t3 ) )
ek_phd15 = np.concatenate( ( ek, ek2, ek3 ) )

################################################

data = np.loadtxt('phd16_ek-1.out')
t = data[:,0]
t = (t*dt)  / (86400*360 )
ek = data[:,5]

data2 = np.loadtxt('phd16_ek-2.out')
t2  = data2[:,0]
ek2 = data2[:,5]
f = np.where(t2 >= 400000)
t2[f], ek2[f] = np.nan, np.nan
t2 = ((t2*dt)  / (86400*360 )) + t[-1]

data3 = np.loadtxt('phd16_ek-3.out')
t3  = data3[:,0]
ek3 = data3[:,5]
f = np.where(t3 >= 400000)
t3[f], ek3[f] = np.nan, np.nan
t3 = ((t3*dt)  / (86400*360 )) + t2[-1]

t_phd16 = np.concatenate( ( t, t2, t3 ) )
ek_phd16 = np.concatenate( ( ek, ek2, ek3 ) )


################################################

data = np.loadtxt('phd18_ek.out')
t = data[:,0]
t_phd18 = (t*dt)  / (86400*360 )
ek_phd18 = data[:,5]

################################################


plt.figure(figsize=(15, 7), facecolor='w')

plt.plot(t_phd15, ek_phd15, 'r', label='Control')
plt.plot(t_phd16, ek_phd16, 'g', label='Flat Bottom')
plt.plot(t_phd18, ek_phd18, 'b', label='BC-only')
plt.grid('on')
plt.axis([0, 10, 0, 0.0015])
plt.legend()
plt.title('Total KE temporal evolution')
plt.xlabel('Years of simulation')
plt.ylabel('J m$^{-2}$')

plt.show()


	

