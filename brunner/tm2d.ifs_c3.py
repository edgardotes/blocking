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
### IFS_4.4-FESOM_5
#xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/z500/daily/IFS_4.4-FESOM_5_z500_day.nc')
### IFS_28-NEMO_25
xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_28-NEMO_25/z500/daily/IFS_28-NEMO_25_z500_day.nc')


### IFS_4.4-FESOM_5
#outpath = '/work/bm1235/b382006/regridding/cycle3/IFS_4.4-FESOM_5/block/ANO_Z500'
### IFS_28-NEMO_25
outpath = '/work/bm1235/b382006/regridding/cycle3/IFS_28-NEMO_25/block/ABS_Z500'

###==subfix 4.4 km  fesom 28 km
#subfix='-ifs-4km'
subfix='-ifs-28km'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'

OUTPATH=outpath+'/'+outfile_flag

###change to yyyy mm dd
xr_in['time']=xr_in.indexes['time'].normalize()
print(xr_in['time'])
### rename zfull
#xr_in=xr_in.rename({'z':'GeopotentialHeight'})



###=== intitate
blk = Blocking()

blk.import_xarray(xr_in)

blk.calculate_gph_from_gp() # calculate geopotential height

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



