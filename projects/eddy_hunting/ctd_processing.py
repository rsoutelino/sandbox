#!/usr/local/epd-7.2-1-rh5-x86_64/bin/python
# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sp
from matplotlib.mlab import griddata
from matplotlib import delaunay
import Image
import glob
from cookb_signalsmooth import smooth
################################################################################

class ReadCTD(object):
	"""docstring for ReadXBT"""
	def __init__(self, filename):
		self.filename = filename
		f = open(filename)
		data = f.readlines()
		data = data[10:]
		self.temp, self.pres, self.cond = [], [], []
		for line in data:
			self.pres.append( float(line.split('\t')[0]) )
			self.temp.append( float(line.split('\t')[1]) )
			self.cond.append( float(line.split('\t')[2]) )
		self.pres = np.array(self.pres, dtype=np.float64)
		self.temp = np.array(self.temp, dtype=np.float64)
		self.cond = np.array(self.cond, dtype=np.float64)

	def find_stations(self, cmin=49, pmin=3, badFlag=200, tips=20):
		"""
		   INPUTS: [integers]
				cmin    : minimum condutivity value to consider water sample 
				pmin    : minimum pressure value to consider water sample
				badFlag : minimum number of air samples to separate stations
				tips    : remove boundary data, related to stabilization of the sensors
		   OUTPUTS: [lists of lists]
		   		returns lists self.PRES, self.TEMP, self.COND as class attributes	

		"""
		count = -1
		bad = []
		self.TEMP, self.PRES, self.COND = [], [], []

		for k in range(self.temp.size):
			if self.cond[k] < cmin:
				# print "Condutivity is ZERO at sample %s" %k
				bad.append(k)
				pass
			else:
				try:
					if self.cond[k-1] < cmin and len(bad) > badFlag and self.pres[k] < pmin:
						count += 1
						self.TEMP.append([])
						self.PRES.append([])
						self.COND.append([])
						# print "\n\nNEW STATION !!! at sample %s\n\n" %k
						# K.append(k)
					else:
						# print "Appending sample %s" %k
						self.TEMP[count].append(self.temp[k])
						self.PRES[count].append(self.pres[k])
						self.COND[count].append(self.cond[k])
						bad = []
				except IndexError: pass

		for k in range(len(self.TEMP)):
			self.PRES[k] = self.PRES[k][tips:-tips] 
			self.TEMP[k] = self.TEMP[k][tips:-tips] 
			self.COND[k] = self.COND[k][tips:-tips]	


	def remove_extremes(self, Prange=(0, 6000), Trange=(0, 40), Crange=(40, 60)):
		"""
		Prange  : (Pmin, Pmax) typical values [tuple]
		Trange  : (Tmin, Tmax) typical values [tuple]
		Crange  : (Cmin, Cmax) typical values [tuple]
		"""
		self.badPres, self.badTemp, self.badCond = [], [], []

		for k in range(len(self.PRES)):
			for l in range(len(self.PRES[k])):
				v1, v2, v3 = self.PRES[k][l], self.TEMP[k][l], self.COND[k][l]
				if (v1 < Prange[0] or v1 > Prange[-1]) or \
				   (v2 < Trange[0] or v2 > Trange[-1]) or \
				   (v3 < Crange[0] or v3 > Crange[-1]):
					self.badPres.append(self.PRES[k][l]); self.PRES[k][l] = None
					self.badTemp.append(self.TEMP[k][l]); self.TEMP[k][l] = None
					self.badCond.append(self.COND[k][l]); self.COND[k][l] = None 

		varlist = [self.PRES, self.TEMP, self.COND]
		for var in varlist:
			for station in var:
				while True:
					try:
						station.remove(None)
					except ValueError:
						break


	def filter_stations(self):
		pass


	def save_raw_stations_in_file(self):
		pass


	def save_filtered_stations_in_file(self):
		pass


################################################################################

ctd = ReadCTD('ctd/V000002.TXT')
print "Finding Stations"
ctd.find_stations()
ctd.remove_extremes()
ctd.save_raw_stations_in_file()
ctd.filter_stations()
ctd.save_filtered_stations_in_file()







			
