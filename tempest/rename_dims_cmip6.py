import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys, os, argparse


fname="cmip6_z500_day.nc"
inpath='/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/daily/'
xr_in=xr.open_dataset(inpath+fname)
print('open..')

### get rid of plev dimmension
##for cycle 2
#gh=xr_in.isel(plev=0,drop=True)

###==rename
print('rename..')
#gh=xr_in.rename({'longitude':'lon',
#                       'latitude':'lat'
#                       })
#

## rename zfull
gh=xr_in.rename({'zg500':'z'})


#
#gh.attrs['units'] = 'm'
#gh.attrs['long_name']= 'Geopotential Height Anomaly'
#
#print('divide g..')

### cts for computing geopotential height
#g = 9.80665  # m s**-2
#geo=gh.z/g

gh.to_netcdf(inpath+'post_'+fname+'')

print("saved")
