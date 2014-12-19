#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seawater.csiro as sw
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import datetime as dt
import glob
import Image

################################################################################

class Paddle(object):
	"""
		INPUT:
			lon:       array or list with longitudes of the track
			lat:       array or list with longitudes of the track
			date:      datetime object with the date of the track
			timedelta: timedelta object with the duration of the track
	"""
	def __init__(self, filename):
		self.filename = filename
		array = np.loadtxt(filename)
		year, month, day = int(filename[:4]), int(filename[4:6]), int(filename[6:8])
		self.date = dt.datetime(year, month, day)
		# self.duration = dt.datetime(0, 0, 0, 0, 0, array[0,1])
		self.track = array[1:,:]
		self.length = self.computeLength()

	def computeLength(self):
		pass


m = Basemap(projection='cyl', llcrnrlat=-23.1, urcrnrlat=-22.9,
             llcrnrlon=-42.1, urcrnrlon=-41.9, lat_ts=0, resolution='f')