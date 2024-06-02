import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys, os, argparse

def cyc_3_ngc3028():
    ### cycle 3 ngc3028
    #fname="ngc3028_z500_day.nc"
    fname="z500_icon_c3_1D*.nc"
    oname="ngc3028_z500_day.nc"
    inpath="/scratch/b/b382006/nextgems/cycle3/ngc3028/z500/daily/"

#    xr_in=xr.open_dataset(inpath+fname,decode_cf=False)
    xr_in=xr.open_mfdataset(inpath+fname)


###** to test
###change to yyyy mm dd
#    xr_in['time']=xr_in.indexes['time'].normalize()
###**

#    xr_in['time'].attrs['units']= 'hours since 2020-01-20'
#    xr_in['lat']=xr_in['lat'][::-1] ## use this because otherwise lats are flipped; only comment for testing

#    gh = xr_in.where(xr_in['zg'] < 9999999.)

## rename zfull
    gh=xr_in.rename({'zg':'z'})

    gh.z.attrs['long_name']= 'Geopotential Height'

#    gh.lat.attrs['standard_name']= 'latitude'
#    gh.lat.attrs['long_name']= 'latitude'
#    gh.lat.attrs['units']= 'degrees_north'
#
#    gh.lon.attrs['standard_name']= 'longitude'
#    gh.lon.attrs['long_name']= 'longitude'
#    gh.lon.attrs['units']= 'degrees_east'

    gh.to_netcdf(inpath+'post_'+oname+'')

def cyc_2_ngc2013():
    ### cycle 2 ngc2013
    fname="ngc2013_z500_day.nc"
    inpath="/scratch/b/b382006/nextgems/cycle2/ngc2013/z500/daily/"

    xr_in=xr.open_dataset(inpath+fname)
    print('open..')

### get rid of plev dimmension
    gh=xr_in.isel(plev=0,drop=True)

## rename zfull
    gh=gh.rename({'zfull':'z'})

    gh.to_netcdf(inpath+'post_'+fname+'')

    print("saved")

def cyc_3_IFS_4_FESOM_5():
    ### cycle 2 ngc2013
    fname="IFS_4.4-FESOM_5_z500_day.nc"
    inpath="/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/z500/daily/"

    xr_in=xr.open_dataset(inpath+fname, decode_cf=False)
    print('open..')

    xr_in['time'].attrs['units']= 'hours since 2020-01-20'

### cts for computing geopotential height
    g = 9.80665  # m s**-2
    gh=xr_in.z/g

    gh.attrs['units'] = 'm'
    gh.attrs['long_name']= 'Geopotential Height'

    gh.to_netcdf(inpath+'post_'+fname+'')

    print("saved")

def cyc_3_IFS_28_NEMO_25():
    ### cycle 3
    fname="IFS_28-NEMO_25_z500_day.nc"
    inpath="/scratch/b/b382006/nextgems/cycle3/IFS_28-NEMO_25/z500/daily/"

    xr_in=xr.open_dataset(inpath+fname, decode_cf=False)
    print('open..')

### get rid of plev dimmension
#    gh=xr_in.isel(plev=0,drop=True)


## rename zfull
#    gh=gh.rename({'zfull':'z'})


    xr_in['time'].attrs['units']= 'hours since 2020-01-20'
#
#    xr_in.z.attrs['units'] = 'm'
#    xr_in.z.attrs['long_name']= 'Geopotential Height'
#
#print('divide g..')

### cts for computing geopotential height
    g = 9.80665  # m s**-2
    gh=xr_in.z/g

    gh.attrs['units'] = 'm'
    gh.attrs['long_name']= 'Geopotential Height'

    gh.to_netcdf(inpath+'post_'+fname+'')

    print("saved")
# ======================================================================================================================================
###### SELECT FUNTION TO RUN ###############

#cyc_2_ngc2013()

cyc_3_ngc3028()

#cyc_3_IFS_4_FESOM_5()

#cyc_3_IFS_28_NEMO_25()
