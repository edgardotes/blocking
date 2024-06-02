import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys, os, argparse


fname="era5_z500_day.nc"

xr_in=xr.open_dataset("/scratch/b/b382006/obs/era5/daily/z500/"+fname+"")
print('open..')



###==rename
gh=xr_in.rename({'longitude':'lon',
                       'latitude':'lat'
                       })

print('rename..')

gh.attrs['units'] = 'm'
gh.attrs['long_name']= 'Geopotential Height Anomaly'

print('divide g..')

### cts for computing geopotential height
g = 9.80665  # m s**-2
geo=gh.z/g

geo.to_netcdf('/scratch/b/b382006/obs/era5/daily/z500/post_'+fname+'')

print("saved")
