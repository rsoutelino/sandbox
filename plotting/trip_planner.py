import datetime as dt
import requests
import matplotlib.pyplot as plt
import numpy as np
from dataIO_tools.read import uds2pandas
from pymo.core.conversions import uv2spddir, k2c



def ms2kts(spd):
    return spd * 1.94384

def uds_request(query, filename):
    # print 'Doing UDS query for {f}: {q}'.format(f=filename, q=query)
    req = requests.get(query)
    with open(filename, 'wb') as f:
        f.write(req.content)

def format_df(df, varmap):
	df = df.rename(columns=varmap)
	df['wsp'], df['wd'] = uv2spddir(df.ugrd10m, df.vgrd10m)
	df.wsp = ms2kts(df.wsp)
	df.gst = ms2kts(df.gst)
	df.tmpsfc = k2c(df.tmpsfc)
	return df

# SETTINGS ##########################################################################
UDS = 'http://metocean:qEwkuwAyyv4iXUEA@uds.metoceanapi.com/uds'
download = False
horizon = 5

sites = dict(
	         tauranga     =  dict(lon=176.1850, lat=-37.6202),
             raglan       =  dict(lon=174.8021, lat=-37.7962),
             newplymouth  =  dict(lon=173.8092, lat=-39.0960),
			 gisborne     =  dict(lon=178.0870, lat=-38.6903)
             )

varmap = {
	       'hs[m]'         :  'hs', 
           'tp[s]'         :  'tp',
		   'dpm[deg]'      :  'dpm',
		   'ugrd10m[m/s]'  :  'ugrd10m',
		   'vgrd10m[m/s]'  :  'vgrd10m',
		   'gst[m/s]'      :  'gst',
		   'tmpsfc[K]'     :  'tmpsfc',
		   'sst[C]'        :  'sst',
		   'sst[]'         :  'sst'
         }	  

plotdict = {
	        'hs': {'max': 4, 'min': 0, 'inc': 0.01, 'size': 30, 'cmap': plt.cm.Blues},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
			# 'hs': {'max': 4, 'min': 0, 'inc': 0.05},
           }

query = "{s}?fmt=txt&var={v}&time={t1},{t2}&xy={x},{y}"
filename_template = '/tmp/{s}_forecast.txt'

#####################################################################################

now = dt.datetime.now()
t1 = dt.datetime(now.year, now.month, now.day)
t2 = t1 + dt.timedelta(days=horizon)

if download:
	for siteID, site in sites.items():
		print "Fetching forecast for {s}".format(s=siteID.upper())
		uds_request(query.format(s=UDS, t1=t1.strftime('%Y%m%d.%H%M%S'), t2=t2.strftime('%Y%m%d.%H%M%S'), 
		                         x=site['lon'], y=site['lat'], v=','.join(varmap.values())), 
								 filename_template.format(s=siteID))

plt.figure(figsize=(11.69, 8.27)) # A4 landscape (fits only 4 sites)
c = 0 

for site in sites.keys():
	c += 1 
	df = format_df(uds2pandas(filename_template.format(s=site)), varmap)
	times = df.index.values
	ax = plt.subplot(2,2,c)

	for idx, time in enumerate(times):
		plot = plotdict['hs']
		hs = np.arange(0, df.hs[idx] - plot['inc'] * 4, plot['inc'])
		df.hs.plot(ax=ax, color='k', linewidth=1)
		plt.scatter([time for v in hs], hs, s=plot['size'], c=hs, marker='s', 
		                                    cmap=plot['cmap'], vmin=plot['min'], 
											vmax=plot['max'], edgecolors=None, alpha=0.5)
		ax.set_ylim([plot['min'], plot['max']])
    
plt.show()



# OLD SIODOC PLOT ####################################################################################

# t = np.arange(0, 24, 0.1)
# m = (np.random.randn(t.size) * 5) + 10
# mraj = m + np.random.rand(t.size)*5 + 1
# m = smooth(m, window_len=15)[:-14]
# mraj = smooth(mraj, window_len=10)[:-9]
# yg = np.arange(0, 25, 0.1)
# xg, yg = np.meshgrid(t, yg)
# direc = ((45*np.pi) / 180) + m*0
# u, v = -1*(m*0+4) , -1*(m*0+4)


# tfill  = np.concatenate( ( np.array([0]), t, np.array([t[-1]]) ) )
# mfill  = np.concatenate( ( np.array([yg.max()]), mraj, np.array([yg.max()]) ) )
# mfill2 = np.concatenate( ( np.array([yg.max()]), m, np.array([yg.max()]) ) )


# fig = plt.figure(facecolor='w', figsize=(14, 10))

# # WIND
# p1 = fig.add_subplot(4,1,1)
# p1.contourf(xg, yg, yg, np.arange(5, 20, 0.1), cmap=plt.cm.hot_r)
# p1.plot(t, m, 'k', linewidth=2)
# p1.fill(tfill, mfill, 'w', edgecolor='w')
# p1.fill(tfill, mfill2, 'w', edgecolor='w', alpha=0.5)
# d = 3
# p1.quiver(t[::d], m[::d]-2, u[::d], v[::d], scale=300, width=0.002, headlength=4.5)
# for k in np.arange(0, t.size, d):
# 	p1.text(t[k], m[k]-1.5, "%0.0f" %m[k], fontsize=8)

# for k in np.arange(0, t.size, d):
# 	p1.text(t[k], mraj[k]+1, "%0.0f" %m[k], color='0.5', fontsize=8)

# p1.axis([0, 10, 0, 25])
# p1.set_axis_off()



# # WAVE
# per = smooth(m, window_len=15)[:-14]

# p2 = fig.add_subplot(4,1,2)
# p2.contourf(xg, yg, yg, np.arange(5, 25, 0.1), cmap=plt.cm.Blues, alpha=0.4)
# p2.plot(t, per, 'r--', linewidth=2)
# p2.plot(t, m*2, 'k', linewidth=2)
# p2.fill(tfill, mfill2*2, 'w', edgecolor='w')

# d = 3
# p2.quiver(t[::d], m[::d]*2-2.3, u[::d]*0, v[::d]*-1, scale=300, width=0.002, headlength=4.5)
# for k in np.arange(0, t.size, d):
# 	p2.text(t[k], m[k]*2+1, "%0.0f" %(m[k]/3), fontsize=8)

# for k in np.arange(0, t.size, d):
# 	p2.text(t[k], per[k]+1, "%0.0f" %per[k], fontsize=8, color='r')

# p2.axis([0, 10, 0, 25])
# p2.set_axis_off()


# # SST
# temp = (np.random.randn(t.size) * 4) + 23
# temp = smooth(temp, window_len=15)[:-14]
# yg = np.arange(12, 28, 0.1)
# xg, yg = np.meshgrid(t, yg)
# tempfill  = np.concatenate( ( np.array([yg.max()]), temp, np.array([yg.max()]) ) )

# p3 = fig.add_subplot(4,1,3)
# p3.contourf(xg, yg, yg, np.arange(12, 28, 0.1), cmap=plt.cm.RdBu_r, alpha=0.8)
# p3.plot(t, temp, 'k', linewidth=2)
# p3.fill(tfill, tempfill, 'w', edgecolor='w')

# d = 3
# for k in np.arange(0, t.size, d):
# 	p3.text(t[k], temp[k]+1, "%0.0f" %temp[k], fontsize=8)

# p3.axis([0, 10, 0, 25])
# p3.set_axis_off()


# # CURRENTS
# p4 = fig.add_subplot(4,1,4)
# p4.contourf(xg, yg, yg, np.arange(12, 28, 0.1), cmap=plt.cm.Greens, alpha=0.8)
# p4.plot(t, temp, 'k', linewidth=2)
# p4.fill(tfill, tempfill, 'w', edgecolor='w')

# d = 3
# for k in np.arange(0, t.size, d):
# 	p4.text(t[k], temp[k]+1, "%0.0f" %temp[k], fontsize=8)

# p4.axis([0, 10, 0, 25])
# p4.set_axis_off()



# plt.show()

