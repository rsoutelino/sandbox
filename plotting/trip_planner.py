import datetime as dt
from collections import OrderedDict
import requests
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dataIO_tools.read import uds2pandas
from pymo.core.conversions import uv2spddir, spddir2uv, k2c



def ms2kts(spd):
    return spd * 1.94384

def uds_request(query, filename):
    # print 'Doing UDS query for {f}: {q}'.format(f=filename, q=query)
    req = requests.get(query)
    with open(filename, 'wb') as f:
        f.write(req.content)

def create_df(df, varmap):
	df = df.rename(columns=varmap)
	df['wsp'], df['wd'] = uv2spddir(df.ugrd10m, df.vgrd10m)
	df['uwave'], df['vwave'] = spddir2uv(df.hs.values * 0 + 1, df.dpm.values)
	df['uwnd'], df['vwnd'] = spddir2uv(df.wsp.values * 0 + 1, df.wd.values)
	df.wsp = ms2kts(df.wsp)
	df.gst = ms2kts(df.gst)
	df.tmpsfc = k2c(df.tmpsfc)
	df = df.rolling(6, win_type='hamming', center=True, closed='both').mean().dropna()
	return df

def create_ax(fig, x, y, w, h, params):
    ax = fig.add_axes([x, y, w, h])
	ax.set_ylim([params['min'], params['max']])
	ax.set_xticks([])
	ax.set_yticks([])
	return ax

def filled_plot(df, var, ax, params):
	times = df.index.values

    for idx, time in enumerate(times):
		values = np.arange(params['min'], df[var][idx] - params['inc'] * 4, params['inc'])
		ax.scatter([time for v in values], values, s=params['size'], c=values, marker='s', 
										   cmap=params['cmap'], vmin=params['min'], 
										   vmax=params['max'], edgecolors=None, alpha=0.5)

def quiver_plot(df, ax, uname, vname, params, qv):
    times = df.index.values
	ax.quiver(times[::['d']], df[var][::['d']] - plot['inc'] * 2, 
	          df[uname][::['d']], df[vname][::['d']], 
			  scale=qv['qsc'], width=qv['qwd'], headlength=qv['hl'])

def annotations(df, ax, var, params, qv):
	times = df.index.values

	for idx in np.arange(0, times.size, qv['d']):
		val = df[var][idx]
	    ax.text(times[idx], val + params['inc'] * 2, '{:0.0f}'.format(val), fontsize=8)

def annotated_scatter(df, ax, var, params, qv):
	times = df.index.values[::qv['d']]
	values = df[var].values[::qv['d']]
	ax.scatter(times, values, s=params['size'], c=values,
	           cmap=params['cmap'], vmin=params['min'], 
			   vmax=params['max'], edgecolors='k', alpha=0.5)

	for idx in np.arange(0, times.size):
		val = values[idx]
	    ax.text(times[idx], val + params['inc'] * 2, '{:0.0f}'.format(val), fontsize=8)

# SETTINGS ##########################################################################

UDS = 'http://metocean:qEwkuwAyyv4iXUEA@uds.metoceanapi.com/uds'
download = False
horizon = 5

OrderedDict([('apple', 4), ('banana', 3), ('orange', 2), ('pear', 1)])

sites = OrderedDict([
                       ('raglan',       dict(lon=174.8021, lat=-37.7962)),	                   
	                   ('tauranga',     dict(lon=176.1850, lat=-37.6202)),
                       ('newplymouth',  dict(lon=173.8092, lat=-39.0960)),
			           ('gisborne',     dict(lon=178.0870, lat=-38.6903))
                    ])

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

plotparams = {
	        'hs':     {'max': 4,  'min': 0,  'inc': 0.05, 'size': 30, 'cmap': plt.cm.Purples},
			'wsp':    {'max': 30, 'min': 0,  'inc': 0.1,  'size': 30, 'cmap': plt.cm.jet},
			'tmpsfc': {'max': 30, 'min': 15, 'inc': 0.05, 'size': 50, 'cmap': plt.cm.jet},
           }

query = "{s}?fmt=txt&var={v}&time={t1},{t2}&xy={x},{y}"
filename_template = '/tmp/{s}_forecast.txt'

# plot settings -------------------

plot_width = 18
inches_per_site = 5
plot_cols = 2.
ax_height = dict(wind=0.2, wave=0.2, weather=0.3)
ax_height_weather = 0.3
xpad = 0.05 
ypad = 0.05 
ax_width = (1 - ypad * 3) / 2

# quiver - magnitude needs to be always 1!
qv = dict(
           d   = 3 # subsampling for quiver
           qsc = 300 # scale
           qwd = 0.002 # width
           hl  = 4.5 # head length
         )
# ---------------------------------

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

# PLOTTING ###########################################################################

nrows = int(ceil(len(sites) / plot_cols))
figsize = (plot_width, nrows * inches_per_site)

fig = plt.figure(figsize=figsize)
c = 0
cur_y_pos = 1 - (ypad + ax_height_wave) # from top to bottom

for site in sites.keys():
	c += 1 
	df = create_df(uds2pandas(filename_template.format(s=site)), varmap)
	times = df.index.values

	if c % 2 == 0:
		x = xpad * 2 + ax_width
	else:
		x = xpad.copy()

	# WIND ----------------------------------------------------
	var = 'wsp'
	params = plotparams[var]
	ax = create_ax(fig, x, y, ax_width, ax_height['wind'], params)
	filled_plot(df, var, ax, params)
    df[var].plot(ax=ax, color='k', linewidth=1)
	quiver_plot(df, ax, 'uwnd', 'vwnd', params, qv)
    annotations(df, ax, params, qv)
	y -= ax_height['wave']

	# WAVE ----------------------------------------------------
    var = 'hs'
	params = plotparams[var]
	ax = create_ax(fig, x, y, ax_width, ax_height['wave'], params)
	filled_plot(df, var, ax, params)
    df[var].plot(ax=ax, color='k', linewidth=1)
	quiver_plot(df, ax, 'uwave', 'vwave', params, qv)
    annotations(df, ax, params, qv)
	y += ax_height['weather']

	# WEATHER -------------------------------------------------
    var = 'tmpsfc'
	params = plotparams[var]
	ax = create_ax(fig, x, y, ax_width, ax_height['weather'], params)
	annotated_scatter(df, var, ax, params)
    df[var].plot(ax=ax, color='k', linewidth=1)
	quiver_plot(df, ax, 'uwave', 'vwave', params, qv)
    annotations(df, ax, params, qv)
	
	y += ax_height['weather']
    











    
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

