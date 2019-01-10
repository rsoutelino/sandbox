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
	df['uwave'], df['vwave'] = df.uwave * -1, df.vwave *- 1
	df['uwnd'], df['vwnd'] = spddir2uv(df.wsp.values * 0 + 1, df.wd.values)
	df.wsp = ms2kts(df.wsp)
	df.gst = ms2kts(df.gst)
	df.tmpsfc = k2c(df.tmpsfc)
	df = df.rolling(6, win_type='hamming', center=True,
	                closed='both').mean().dropna()
	return df


def create_ax(fig, x, y, w, h, params):
    ax = fig.add_axes([x, y, w, h])
    ax.set_ylim([params['min'], params['max']])
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


def filled_plot(df, var, ax, params):
	# import pdb; pdb.set_trace()
	times = df.index
	ax.set_xlim([df.index[0], df.index[-1]])
	for idx, time in enumerate(times):
		values = np.arange(params['min'], df[var][idx] - params['inc'] * 4, params['inc'])
        ax.scatter([time for v in values], values, s=params['size'], c=values,
		           marker='s', cmap=params['cmap'], 
				   vmin=params['min'], vmax=params['max'], 
				   edgecolors=None, alpha=0.5)


def quiver_plot(df, var, ax, uname, vname, params, qv):
	times = df.index.values[::qv['d']]
	y = df[var][::qv['d']]
	u = df[uname][::qv['d']]
	v = df[vname][::qv['d']]
	ax.set_xlim([df.index[0], df.index[-1]])
	ax.quiver(times, y - params['inc'] * 10, u, v, 
			  scale=qv['qsc'], width=qv['qwd'], 
              headlength=qv['hl'], pivot='middle')


def annotations(df, var, ax, params, qv, fmt, **kwargs):
	times = df.index.values

	for idx in np.arange(0, times.size, qv['d']):
		val = df[var][idx]
		ax.text(times[idx], val + params['inc'] * 5, fmt.format(val), fontsize=8, **kwargs)


def annotated_scatter(df, var, ax, params, qv):
	times = df.index.values[::qv['d']]
	values = df[var].values[::qv['d']]
	ax.scatter(times, values, s=params['size'], c=values,
	           cmap=params['cmap'], vmin=params['min'],
			   vmax=params['max'], edgecolors='k', alpha=0.5)

	ax.set_xlim([df.index[0], df.index[-1]])

	for idx in np.arange(0, times.size):
		val = values[idx]
		ax.text(times[idx], val + params['inc'] * 2, '{:0.0f}'.format(val), 
                fontsize=6, horizontalalignment='center', verticalalignment='center')


# SETTINGS ##########################################################################

toff = '13H'
UDS = 'http://metocean:qEwkuwAyyv4iXUEA@uds.metoceanapi.com/uds'
download = False
horizon = 5

OrderedDict([('apple', 4), ('banana', 3), ('orange', 2), ('pear', 1)])

sites = OrderedDict([
                       ('raglan',       dict(lon=174.8021, lat=-37.7962)),	                   
	                   ('tauranga',     dict(lon=176.1850, lat=-37.6202)),
                       ('newplymouth',  dict(lon=173.8092, lat=-39.0960)),
			           ('gisborne',     dict(lon=178.0870, lat=-38.6903)),
					   ('ahipara',      dict(lon=173.1088, lat=-35.1289)),
					   ('northland',    dict(lon=174.0976, lat=-35.2232)),
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
			'wsp':    {'max': 30, 'min': 0,  'inc': 0.4,  'size': 30, 'cmap': plt.cm.jet},
			'tmpsfc': {'max': 30, 'min': 15, 'inc': 0.05, 'size': 120, 'cmap': plt.cm.jet},
			'tp':     {'max': 20, 'min': 0,  'inc': 0.05, 'size': None, 'cmap': None},
           }

query = "{s}?fmt=txt&var={v}&time={t1},{t2}&xy={x},{y}"
filename_template = '/tmp/{s}_forecast.txt'

# plot settings -------------------

plot_width = 18
inches_per_site = 5
plot_cols = 2.
ax_height = dict(wind=0.1, wave=0.1, weather=0.1)
ax_height_total = ax_height['wind'] + ax_height['wave'] + ax_height['weather']
xpad = 0.02 
ypad = 0.02 
ax_width = (1 - ypad * 3) / 2

# quiver - magnitude needs to be always 1!
qv = dict(
           d   = 6, # subsampling for quiver
           qsc = 60, # scale
           qwd = 0.002, # width
           hl  = 4.5, # head length
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

nrows = int(np.ceil(len(sites) / plot_cols))
figsize = (plot_width, nrows * inches_per_site)

fig = plt.figure(figsize=figsize)
c = 0
ypos = 1 - (ypad + ax_height['wind']) # from top to bottom

# import pdb; pdb.set_trace()
for site in sites.keys():
	c += 1 
	df = create_df(uds2pandas(filename_template.format(s=site)), varmap)
	df.index += pd.to_timedelta(toff) # timezone
	times = df.index.values

	if c % 2 == 0:
		xpos = xpad * 2 + ax_width
	else:
		xpos = xpad

	# WIND ----------------------------------------------------
	var = 'wsp'
	params = plotparams[var]
	ax = create_ax(fig, xpos, ypos, ax_width, ax_height['wind'], params)
	ax.set_title(site.capitalize(), fontweight='bold', fontsize=8)
	# filled_plot(df, var, ax, params)
	df[var].plot(ax=ax, color='k')
	quiver_plot(df, var, ax, 'uwnd', 'vwnd', params, qv)
	annotations(df, var, ax, params, qv, '{:0.0f}')
	ypos -= ax_height['wave']

	# WAVE ----------------------------------------------------
	var = 'hs'
	params = plotparams[var]
	ax1 = create_ax(fig, xpos, ypos, ax_width, ax_height['wave'], params)
	# filled_plot(df, var, ax, params)
	df[var].plot(ax=ax1, color='k')
	quiver_plot(df, var, ax1, 'uwave', 'vwave', params, qv)
	annotations(df, var, ax1, params, qv, '{:0.1f}')
    # ---------------------------------------------------------
	var = 'tp'
	params = plotparams[var]
	ax2 = create_ax(fig, xpos, ypos, ax_width, ax_height['wave'], params)
	df[var].plot(ax=ax2, color='r', dashes=[3, 3])
	annotations(df, var, ax2, params, qv, '{:0.0f}', color='r')

	ypos -= ax_height['weather']

	# WEATHER -------------------------------------------------
	var = 'tmpsfc'
	params = plotparams[var]
	ax = create_ax(fig, xpos, ypos, ax_width, ax_height['weather'], params)
	annotated_scatter(df, var, ax, params, qv)

    # spotting the ycoord of the next site
	if c % 2 == 0:
		ypos -= (ypad + ax_height['wind'])
	else:
		ypos += (ax_height_total - ax_height['wind'])


plt.show()
