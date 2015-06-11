#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import time


imp = 'pg'
rootdir = os.path.join("/ops/forecast/roms", imp)

filelist = glob.glob("%s/*.out")

plt.figure()


dt = 300

os.system('sh roms_ek-3.sh')
data = np.loadtxt('phd16_ek-3.out')
t = data[:,0]
t = (t*dt)  / (86400*360 )
ek = data[:,5]
plt.plot(t, ek, 'b')
plt.grid('on')
plt.title('Total KE temporal evolution')
plt.xlabel('Years of simulation')
plt.ylabel('J m$^{-2}$')
plt.show()


	

