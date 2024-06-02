import matplotlib.pylab as plt
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import cm
import cartopy.crs as ccrs

f = xr.open_dataset('/Users/edgardolorestesillos/Documents/NEXTGEMS/tools/blocking-master/data/geop_1degree_bf_3D.nc')

#f.Blocking.plot()

acc=f.Blocking.sum(dim='time')

###=== plotting
ax = plt.axes(projection=ccrs.NorthPolarStereo())
acc.plot.contourf(levels = np.linspace(0,30,10),cmap = cm.hot_r,ax=ax, transform=ccrs.PlateCarree())

ax.coastlines()
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())

plt.show()
