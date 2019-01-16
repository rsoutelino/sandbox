import datetime as dt
from collections import OrderedDict
import requests
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dataIO_tools.read import uds2pandas
from pymo.core.conversions import uv2spddir, spddir2uv, k2c
from astral import Astral
import uuid

def ms2kts(spd):
    return spd * 1.94384


def uds_request(query, filename):
    # print 'Doing UDS query for {f}: {q}'.format(f=filename, q=query)
    req = requests.get(query)
    with open(filename, 'wb') as f:
        f.write(req.content)


def create_df(df, varmap, toff):
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
    df.index += pd.to_timedelta(toff) # timezone
    return df


def create_ax(df, fig, x, y, w, h, params):
    ax = fig.add_axes([x, y, w, h], label=str(uuid.uuid1()))
    ax.set_ylim([params['min'], params['max']])
    ax.set_xlim([df.index[0], df.index[-1]])
    ax.axis('off')
    return ax


def plot_nightshade(df, ax, **kwargs):
    a = Astral()
    a.solar_depression = 'civil'
    city = a['wellington']
    ymin, ymax = ax.get_ylim()
    
    for day in pd.date_range(df.index[0].date(), df.index[-1].date()):
        sun1 = city.sun(date=day - dt.timedelta(days=1))
        sun2 = city.sun(date=day, local=True)
        sunset = sun1['sunset'].replace(tzinfo=None)
        sunrise = sun2['sunrise'].replace(tzinfo=None)
        night = pd.DataFrame(index=[sunset, sunrise], 
		             data=dict(shade=[ymax, ymax]))
        night.shade.plot(kind='area', ax=ax, color='0.9', alpha=0.5, **kwargs)


def filled_plot(df, var, ax, params, **kwargs):
    times = df.index
    ax.set_xlim([df.index[0], df.index[-1]])

    for idx, time in enumerate(times):
        values = np.arange(params['min'], df[var][idx] - params['inc'] * 4, params['inc'])
        ax.scatter([time for v in values], values, s=params['size'], c=values,
		    marker='s', cmap=params['cmap'],
                    vmin=params['min'], vmax=params['max'], 
                    edgecolors=None, alpha=0.5, **kwargs)


def quiver_plot(df, var, ax, uname, vname, params, qv): 
    times = df.index.values[::qv['d']]
    y = df[var][::qv['d']]
    u = df[uname][::qv['d']]
    v = df[vname][::qv['d']]
    ax.set_xlim([df.index[0], df.index[-1]])
    ax.quiver(times, y - params['inc'] * 10, u, v, 
              scale=qv['qsc'], width=qv['qwd'], 
              headlength=qv['hl'], pivot='middle', zorder=2)


def annotations(df, var, ax, params, qv, fmt, skip=None, **kwargs):
    times = df.index.values
    if not skip: 
        skip = qv['d']

    for idx in np.arange(0, times.size, skip):
        val = df[var][idx]
        ax.text(times[idx], val + params['inc'] * 5, fmt.format(val), fontsize=6, **kwargs)
    
    ax.set_xlim([df.index[0], df.index[-1]])


def annotated_scatter(df, var, ax, params, qv, skip=None, **kwargs):
    if not skip:
        skip = qv['d']

    times = df.index.values[::skip]
    values = df[var].values[::skip]
    ax.scatter(times, values, s=params['size'], c=values,
               cmap=params['cmap'], vmin=params['min'],
               vmax=params['max'], edgecolors='k', **kwargs)

    ax.set_xlim([df.index[0], df.index[-1]])

    for idx in np.arange(0, times.size):
        val = values[idx]
        ax.text(times[idx], val + params['inc'] * 2, '{:0.0f}'.format(val), 
                fontsize=6, horizontalalignment='center', 
                verticalalignment='center', zorder=6)


def plot_prate(df, var, ax, params, **kwargs):
    df = df.resample('6H').mean()
    df['precip'] = df.apratesfc * 3600
    df.precip.plot(ax=ax, kind='bar', color='b', alpha=0.5)


def plot_time(df, ax, now, horizon):
    y = ax.get_ylim()[-1]
    t1 = dt.datetime(now.year, now.month, now.day, 12)
    t2 = t1 + dt.timedelta(days=horizon - 1)
    times = pd.date_range(t1, t2)
    
    for time in times:
        txt = "{w} {d}.{m}".format(w=time.strftime('%a'), d=time.day, m=time.month)
        ax.text(time, y - 1.2, txt, fontsize=8, color='k',
                horizontalalignment='center', 
                verticalalignment='center')

    times = pd.date_range(t1, t2, freq='6H')

    for time in times:
        txt = "{h}h".format(h=time.hour)
        ax.text(time, y - 4, txt, fontsize=6, color='0.5', 
                horizontalalignment='center', 
                verticalalignment='center')

    ax.set_xlim(df.index[0], df.index[-1])


# SETTINGS ##########################################################################

toff = '13H'
UDS = 'http://metocean:qEwkuwAyyv4iXUEA@uds.metoceanapi.com/uds'
download = False
horizon = 7

sites = OrderedDict([
                      ('raglan',       dict(lon=174.8021, lat=-37.7962, lon2=174.9015, lat2=-37.8035)),	                   
	                  ('tauranga',     dict(lon=176.1850, lat=-37.6202, lon2=176.1944, lat2=-37.6969)),
                    #    ('newplymouth',  dict(lon=173.8092, lat=-39.0960, lon2=174.1173, lat2=-39.0669)),
                      ('whangamata',  dict(lon=175.9072, lat=-37.2060, lon2=175.8636, lat2=-37.2149)),
		              ('gisborne',     dict(lon=178.0870, lat=-38.6903, lon2=178.0149, lat2=-38.6394)),
					#    ('ahipara',      dict(lon=173.1088, lat=-35.1289, lon2=173.2132, lat2=-35.1137)),
					#    ('keri-keri',    dict(lon=174.0976, lat=-35.2232, lon2=174.0543, lat2=-35.2724)),
                    ])

varmap = {
	       'hs[m]'                :  'hs', 
           'tp[s]'                :  'tp',
		   'dpm[deg]'             :  'dpm',
		   'ugrd10m[m/s]'         :  'ugrd10m',
		   'vgrd10m[m/s]'         :  'vgrd10m',
		   'gst[m/s]'             :  'gst',
		   'tmpsfc[K]'            :  'tmpsfc',
		   'sst[C]'               :  'sst',
		   'sst[]'                :  'sst',
		   'apratesfc[kg/m^2/s]'  :  'apratesfc',
		   'tcdcclm[%]'           :  'tcdcclm',
         }	  

plotparams = {
	        'hs':        {'max': 4,   'min': 0.2,  'inc': 0.03, 'size': 30, 'cmap': plt.cm.Blues},
			'wsp':       {'max': 30,  'min': 5,  'inc': 0.3,  'size': 30, 'cmap': plt.cm.YlOrBr},
			'tmpsfc':    {'max': 35,  'min': 5, 'inc': 0.05, 'size': 120, 'cmap': plt.cm.Spectral_r},
			'tp':        {'max': 23,  'min': 0,  'inc': 0.1, 'size': None, 'cmap': None},
			'apratesfc': {'max': 10,  'min': 0,  'inc': 0.05, 'size': None, 'cmap': None},
			'tcdcclm':   {'max': 100, 'min': 0,  'inc': None, 'size': None, 'cmap': None},			
			'et':        {'max': 5,   'min': 0,  'inc': None, 'size': None, 'cmap': None},			
			'sst':       {'max': 25,  'min': 0, 'inc': 0.05, 'size': 120, 'cmap': plt.cm.Spectral_r},			
           }

query = "{s}?fmt=txt&var={v}&time={t1},{t2}&xy={x},{y}"
tide_query = "{s}?fmt=txt&var=et&time={t1},{t2}&xy={x},{y}&type=tide&datum=LAT&dt=0.1"
filename_template = '/tmp/{s}_forecast.txt'

# plot settings -------------------

plot_width = 18
inches_per_site = 5
plot_cols = 2.
ax_height = dict(wind=0.15, wave=0.15, weather=0.15)
ax_height_total = ax_height['wind'] + ax_height['wave'] + ax_height['weather']
xpad = 0.04 
ypad = 0.04 
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
t1 -= pd.to_timedelta('16H') # because we loose 3H in the rolling mean
t2 -= pd.to_timedelta('10H')

if download:
    for siteID, site in sites.items():
        print "Fetching forecast for {s}".format(s=siteID.upper())
        uds_request(query.format(s=UDS, t1=t1.strftime('%Y%m%d.%H%M%S'), t2=t2.strftime('%Y%m%d.%H%M%S'), 
                                 x=site['lon'], y=site['lat'], v=','.join(varmap.values())), 
								 filename_template.format(s=siteID))
        uds_request(query.format(s=UDS, t1=t1.strftime('%Y%m%d.%H%M%S'), t2=t2.strftime('%Y%m%d.%H%M%S'),
                                 x=site['lon2'], y=site['lat2'], v=','.join(varmap.values())),
                                 filename_template.format(s=siteID + '_land'))
        uds_request(tide_query.format(s=UDS, t1=t1.strftime('%Y%m%d.%H%M%S'), t2=t2.strftime('%Y%m%d.%H%M%S'),
                                 x=site['lon'], y=site['lat'], v=','.join(varmap.values())),
                                 filename_template.format(s=siteID + '_tide'))

# PLOTTING ###########################################################################

nrows = int(np.ceil(len(sites) / plot_cols))
figsize = (plot_width, nrows * inches_per_site)

fig = plt.figure(figsize=figsize)
c = 0
ypos = 1 - (ypad + ax_height['wind']) # from top to bottom

# import pdb; pdb.set_trace()
for site in sites.keys():
    c += 1 
    df = create_df(uds2pandas(filename_template.format(s=site)), varmap, toff)
    df2 = create_df(uds2pandas(filename_template.format(s=site + '_land')), varmap, toff)
    tide = uds2pandas(filename_template.format(s=site + '_tide'))
    tide.index += pd.to_timedelta(toff)
    tide = tide.rename(columns={'et[m]': 'et'})
    times = df.index.values

    if c % 2 == 0:
        xpos = xpad * 2 + ax_width
    else:
        xpos = xpad

    # WIND ----------------------------------------------------
    var = 'wsp'
    params = plotparams[var]
    ax = create_ax(df, fig, xpos, ypos, ax_width, ax_height['wind'], params)
    plot_time(df, ax, now, horizon)
    tit = ax.set_title(site.capitalize(), fontweight='bold', fontsize=10)
    tit.set_position([0.5, 1])
    plot_nightshade(df, ax)
    filled_plot(df, var, ax, params)
    df[var].plot(ax=ax, color='k')
    quiver_plot(df, var, ax, 'uwnd', 'vwnd', params, qv)
    annotations(df, var, ax, params, qv, '{:0.0f}')
    ypos -= ax_height['wave']

    # WAVE ----------------------------------------------------
    var = 'hs'
    params = plotparams[var]
    ax1 = create_ax(df, fig, xpos, ypos, ax_width, ax_height['wave'], params)
    plot_nightshade(df, ax1)
    filled_plot(df, var, ax1, params, zorder=2)
    df[var].plot(ax=ax1, color='k', zorder=2)
    quiver_plot(df, var, ax1, 'uwave', 'vwave', params, qv)
    annotations(df, var, ax1, params, qv, '{:0.1f}')
    # ---------------------------------------------------------
    var = 'tp'
    params = plotparams[var]
    ax = create_ax(df, fig, xpos, ypos, ax_width, ax_height['wave'], params)
    df[var].plot(ax=ax, color='r', dashes=[3, 3], zorder=3)
    annotations(df, var, ax, params, qv, '{:0.0f}', color='r', skip=12)
    # ---------------------------------------------------------
    var = 'sst'
    params = plotparams[var]
    ax = create_ax(df, fig, xpos, ypos, ax_width, ax_height['wave'], params)
    df[var].plot(ax=ax, color='0.6', zorder=4)
    annotated_scatter(df, var, ax, params, qv, zorder=5, skip=24)

    ypos -= ax_height['weather']

    # WEATHER -------------------------------------------------
    var = 'tcdcclm'
    params = plotparams[var]
    ax = create_ax(df2, fig, xpos, ypos, ax_width, ax_height['weather'], params)
    plot_nightshade(df2, ax, zorder=1)
    df2[var].plot(ax=ax, kind='area', color='0.7', alpha=0.5, zorder=2)
    # ---------------------------------------------------------
    var = 'et'
    params = plotparams[var]
    ax = create_ax(tide, fig, xpos, ypos, ax_width, ax_height['weather'], params)
    tide.et.plot(ax=ax, kind='area', color='g', alpha=0.2, zorder=3)
    ax.set_xlim([df.index[0], df.index[-1]])
    # ---------------------------------------------------------
    var = 'tmpsfc'
    params = plotparams[var]
    ax = create_ax(df2, fig, xpos, ypos, ax_width, ax_height['weather'], params)
    df2[var].plot(ax=ax, color='k', zorder=2)
    annotated_scatter(df2, var, ax, params, qv, zorder=4)
    # ---------------------------------------------------------
    var = 'apratesfc'
    params = plotparams[var]
    ax = create_ax(df2, fig, xpos, ypos, ax_width, ax_height['weather'], params)
    plot_prate(df2, var, ax, params, zorder=5)

    # spotting the ycoord of the next site
    if c % 2 == 0:
        ypos -= (ypad + ax_height['wind'])
    else:
        ypos += (ax_height_total - ax_height['wind'])

plt.savefig('/tmp/trip_planner.png')
plt.show()
