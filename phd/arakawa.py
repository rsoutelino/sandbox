import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

grid = nc.Dataset('phd15_grd.nc')

fig = plt.figure(facecolor='w')

plt.plot(grid.variables['lon_rho'][:6,:6], grid.variables['lat_rho'][:6,:6], '0.9' )
for k in range(0,6):
    plt.plot(grid.variables['lon_rho'][k,:6], grid.variables['lat_rho'][k,:6], '0.9' )

plt.plot(grid.variables['lon_rho'][:6,:6], grid.variables['lat_rho'][:6,:6],'b*', markersize=10, label='rho points')
plt.plot(grid.variables['lon_u'][:6,:5], grid.variables['lat_u'][:6,:5],'ro', markersize=5, label='u points')
plt.plot(grid.variables['lon_v'][:5,:6], grid.variables['lat_v'][:5,:6],'go', markersize=5, label='v points')

plt.text(-40.736, -22.772, 'Rho-points', color='blue', fontsize=10, fontweight='bold')
plt.text(-40.736, -22.794, 'U-points', color='red', fontsize=10, fontweight='bold')
plt.text(-40.736, -22.818, 'V-points', color='green', fontsize=10, fontweight='bold')

plt.axis('equal')
plt.axis('off')
plt.title('ROMS Grid set-up')
plt.show()