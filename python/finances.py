import pylab as pl
import numpy as np
import matplotlib as plt
import datetime as dt

### finances.py


### MARCH: settle down month ###############################################

# EQUIPS: Clothes, shoes, jackets, GPS, Cell phones
eq03 = 473.77

# FIXED: house things, CAPES stipends, AAA-YMCA-join fees,
#  phone refil card, car insurance initial fee, car registration and inspection  
fix03 = 837		
	
### APRIL: apt mouting month ###############################################

# EQUIPS: clothes, iPod docking station
eq04 = 156

# FIXED: house things, furniture, etc
fix04 = 2431

#############################################################################
### MONTHLY BUGDET EVALUATION ###############################################
#############################################################################
VAR = np.loadtxt('finances.txt')

years  = np.array( VAR[:,0], dtype=int )
months = np.array( VAR[:,1], dtype=int )
days   = np.array( 1 + 0*years, dtype=int )
fixed  = VAR[:,2]
phone  = VAR[:,3]
util   = VAR[:,4]
groc   = VAR[:,5]
fuel   = VAR[:,6]
cmain  = VAR[:,7]
fun    = VAR[:,8]
equip  = VAR[:,9]

dates = []

for year, month, day in zip(years, months, days):
	dates.append( dt.datetime(year, month, day) )

dates2 = plt.dates.date2num(dates) # matlab-like dates

base   = fixed + util + phone + groc + fuel
extra  = cmain + fun + equip
total  = base + extra

fmt = plt.ticker.Formatter()

#### GRAPHICS
fig1 = pl.figure(1,figsize=(15,9),facecolor='w')
p1 = pl.subplot(311)
p1.plot_date(dates2,base,'k-*',linewidth=2,label='Base')
p1.plot_date(dates2,extra,'r-*',linewidth=2,label='Extra')
p1.plot_date(dates2,total,'g-',linewidth=2,label='TOTAL')
p1.legend()
p1.xaxis.set_major_formatter(plt.dates.DateFormatter('%b'))
p1.grid(b=True)
p1.axis([dt.datetime(2010,2,15), dt.datetime(2011,9,15), -100, 3500])

p2 = pl.subplot(312)
p2.plot_date(dates2,phone,'k--*',linewidth=1,label='Phone')
p2.plot_date(dates2,util,'-*',color=('0.5'),linewidth=2,label='Utilities')
p2.plot_date(dates2,groc,'b-*',linewidth=2,label='Groceries')
p2.plot_date(dates2,fuel,'r-*',linewidth=2,label='Fuel')
p2.legend()
p2.xaxis.set_major_formatter(plt.dates.DateFormatter('%b'))
p2.grid(b=True)
p2.axis([dt.datetime(2010,2,15), dt.datetime(2011,9,15), -10, 500])

p3 = pl.subplot(313)
p3.plot_date(dates2,cmain,'r-*',linewidth=2,label='Car Maint')
p3.plot_date(dates2,fun,'b-*',linewidth=2,label='Fun')
p3.plot_date(dates2,equip,'k-*',color=('0.5'),linewidth=2,label='Equip')
p3.legend()
p3.xaxis.set_major_formatter(plt.dates.DateFormatter('%b'))
p3.grid(b=True)
p3.axis([dt.datetime(2010,2,15), dt.datetime(2011,9,15), -10, 1000])



p3.text(dt.datetime(2010,7,1),550,'Breaks',color=('0.5'))
p3.text(dt.datetime(2010,5,15),400,'Mattress',color=('0.5'))
p3.text(dt.datetime(2010,11,1),800,'PR Trip',color=('0.5'))
p3.text(dt.datetime(2011,1,15),800,'Miami',color=('0.5'))
p3.text(dt.datetime(2011,3,1),800,'Oregon',color=('0.5'))


pl.show()















