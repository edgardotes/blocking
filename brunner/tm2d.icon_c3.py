import sys
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm

sys.path.append("/home/b/b382006/tools/Blocking/brunner/blocking/")
from BlockingDetection import Blocking

### blocking method based on the original Code of Lukas Brunner
#xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/ngc3028/z500/daily/ngc3028_z500_day.nc')
xr_in=xr.open_mfdataset('/scratch/b/b382006/nextgems/cycle3/ngc3028/z500/daily/z500_icon_c3_1D_*nc')

###==ifs
outpath = '/work/bm1235/b382006/regridding/cycle3/ngc3028/block/ABS_Z500/'

###==subfix ngc2013
subfix='-icon-ngc3028'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'

OUTPATH=outpath+'/'+outfile_flag

### rename zfull
xr_in=xr_in.rename({'zg':'GeopotentialHeight'})

###change to yyyy mm dd
#xr_in['time']=xr_in.indexes['time'].normalize()

###Drop dimension time_bnds
#xr_in=xr_in.drop_dims('bnds')

###=== intitate
blk = Blocking()

blk.import_xarray(xr_in)

#blk.calculate_gph_from_gp() # calculate geopotential height

###
blk.set_up(time_name='time', longitude_name='lon', latitude_name='lat')

blk.calculate_gradients(delta_degree=15)

blk.calculate_ib(
    gradient_equator_below=0,
    gradient_pole_below=-10,
    gradient_equator2_above=5)

#
blk.calculate_eib(min_extent_degree=15)

#
blk.calculate_blocking(
    stationary_pm_days=2,
    longitude_pm_degree=8,
    latitude_pm_degree=2)

# save to disk
blk.save(OUTPATH, 'Blocking')



