
import matplotlib.pyplot as plt
import numpy as np
from cookb_signalsmooth  import smooth
from intdir2uv import *




t = np.arange(0, 24, 0.1)
m = (np.random.randn(t.size) * 5) + 10
mraj = m + np.random.rand(t.size)*5 + 1
m = smooth(m, window_len=15)
mraj = smooth(mraj, window_len=10)
yg = np.arange(0, 25, 0.1)
xg, yg = np.meshgrid(t, yg)
direc = ((45*np.pi) / 180) + m*0
u, v = -1*(m*0+4) , -1*(m*0+4)


tfill  = np.concatenate( ( np.array([0]), t, np.array([t[-1]]) ) )
mfill  = np.concatenate( ( np.array([yg.max()]), mraj, np.array([yg.max()]) ) )
mfill2 = np.concatenate( ( np.array([yg.max()]), m, np.array([yg.max()]) ) )


fig = plt.figure(facecolor='w', figsize=(14, 10))

# WIND
p1 = fig.add_subplot(4,1,1)
p1.contourf(xg, yg, yg, np.arange(5, 20, 0.1), cmap=plt.cm.hot_r)
p1.plot(t, m, 'k', linewidth=2)
p1.fill(tfill, mfill, 'w', edgecolor='w')
p1.fill(tfill, mfill2, 'w', edgecolor='w', alpha=0.5)
d = 3
p1.quiver(t[::d], m[::d]-2, u[::d], v[::d], scale=300, width=0.002, headlength=4.5)
for k in np.arange(0, t.size, d):
	p1.text(t[k], m[k]-1.5, "%0.0f" %m[k], fontsize=8)

for k in np.arange(0, t.size, d):
	p1.text(t[k], mraj[k]+1, "%0.0f" %m[k], color='0.5', fontsize=8)

p1.axis([0, 10, 0, 25])
p1.set_axis_off()



# WAVE
per = smooth(m, window_len=15)

p2 = fig.add_subplot(4,1,2)
p2.contourf(xg, yg, yg, np.arange(5, 25, 0.1), cmap=plt.cm.Blues, alpha=0.4)
p2.plot(t, per, 'r--', linewidth=2)
p2.plot(t, m*2, 'k', linewidth=2)
p2.fill(tfill, mfill2*2, 'w', edgecolor='w')

d = 3
p2.quiver(t[::d], m[::d]*2-2.3, u[::d]*0, v[::d]*-1, scale=300, width=0.002, headlength=4.5)
for k in np.arange(0, t.size, d):
	p2.text(t[k], m[k]*2+1, "%0.0f" %(m[k]/3), fontsize=8)

for k in np.arange(0, t.size, d):
	p2.text(t[k], per[k]+1, "%0.0f" %per[k], fontsize=8, color='r')

p2.axis([0, 10, 0, 25])
p2.set_axis_off()


# SST
temp = (np.random.randn(t.size) * 4) + 23
temp = smooth(temp, window_len=15)
yg = np.arange(12, 28, 0.1)
xg, yg = np.meshgrid(t, yg)
tempfill  = np.concatenate( ( np.array([yg.max()]), temp, np.array([yg.max()]) ) )

p3 = fig.add_subplot(4,1,3)
p3.contourf(xg, yg, yg, np.arange(12, 28, 0.1), cmap=plt.cm.RdBu_r, alpha=0.8)
p3.plot(t, temp, 'k', linewidth=2)
p3.fill(tfill, tempfill, 'w', edgecolor='w')

d = 3
for k in np.arange(0, t.size, d):
	p3.text(t[k], temp[k]+1, "%0.0f" %temp[k], fontsize=8)

p3.axis([0, 10, 0, 25])
p3.set_axis_off()


# CURRENTS
p4 = fig.add_subplot(4,1,4)
p4.contourf(xg, yg, yg, np.arange(12, 28, 0.1), cmap=plt.cm.Greens, alpha=0.8)
p4.plot(t, temp, 'k', linewidth=2)
p4.fill(tfill, tempfill, 'w', edgecolor='w')

d = 3
for k in np.arange(0, t.size, d):
	p4.text(t[k], temp[k]+1, "%0.0f" %temp[k], fontsize=8)

p4.axis([0, 10, 0, 25])
p4.set_axis_off()



plt.show()

