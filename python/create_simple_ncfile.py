#!/usr/bin/env python

import datetime as dt

now = dt.datetime.now()

ncfile = nc.Dataset('exemplo.nc', mode='w', clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('x', size=150)
ncfile.createDimension('y', size=150)
ncfile.createDimension('z', size=121)

# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', 'Interpolated T/S Climatology Fields')
setattr(ncfile, 'date', str(now))

# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('lon', 'd', dimensions=('x', 'y'))
setattr(ncfile.variables['lon'], 'long_name', 'longitude')
setattr(ncfile.variables['lon'], 'units', 'degree')
ncfile.variables['lon'][:]  = lon

# ---------------------------------------------------------------------------
ncfile.createVariable('lat', 'd', dimensions=('x', 'y'))
setattr(ncfile.variables['lat'], 'long_name', 'latitude')
setattr(ncfile.variables['lat'], 'units', 'degree')
ncfile.variables['lat'][:]  = lat

# ---------------------------------------------------------------------------
ncfile.createVariable('depth', 'd', dimensions=('z'))
setattr(ncfile.variables['depth'], 'long_name', 'depth')
setattr(ncfile.variables['depth'], 'units', 'meters')
ncfile.variables['depth'][:]  = z

# ---------------------------------------------------------------------------
ncfile.createVariable('temp', 'd', dimensions=('x', 'y', 'z'))
setattr(ncfile.variables['temp'], 'long_name', 'temperature')
setattr(ncfile.variables['temp'], 'units', 'degree celcius')
ncfile.variables['temp'][:]  = temp

# ---------------------------------------------------------------------------
ncfile.createVariable('salt', 'd', dimensions=('x', 'y', 'z'))
setattr(ncfile.variables['salt'], 'long_name', 'salinity')
setattr(ncfile.variables['salt'], 'units', 'none')
ncfile.variables['salt'][:]  = temp

ncfile.sync()
