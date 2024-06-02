import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys
sys.path.append("/home/b/b382006/tools/Blocking/brunner/blocking/")
from BlockingDetection import Blocking

### blocking method based on the original Code of Lukas Brunner
xr_in=xr.open_mfdataset('/work/bm1235/b382006/regridding/ifs/tco1279-orca025_9km_to_1degree/Y*/ICMUAhr2n+0*')

###==ngc2013
outpath = '/work/bm1235/b382006/regridding/ifs/tco1279-orca025_9km_to_1degree/block/ABS_Z500'

###==subfix ngc2013
subfix='-ifs'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'

#Display time
#xr_in.time

### Take geopotential at 500 hPa
in_z_500=xr_in.z.sel(plev=50000.0)

###daily mean
daily_z=in_z_500.resample(time='1D').mean(dim='time',keep_attrs=True)

###select months
###pd_z = daily_z.sel(time=daily_z.time.dt.month.isin([1, 2, 12]))
###dataset_z=xr.DataArray.to_dataset(pd_z)
dataset_z=xr.DataArray.to_dataset(daily_z)


### configure attributes 
dataset_z.time.attrs['standard_name'] = 'time'
dataset_z.time.attrs['axis'] = 'T'
dataset_z.time

### Outpath line
#OUTPATH = dataset_z.replace(
#    'geopotential', 'blocking').replace('.nc', '_bf_3D.nc')
OUTPATH=outpath+'/'+outfile_flag

###=== intitate
blk = Blocking()
blk.import_xarray(dataset_z)

###=== pre-processing
#blk.get_time_subset(months='DJF')
#blk.calculate_daily_mean()
blk.calculate_gph_from_gp() # calculate geopotential height
#blk.ds.time

blk.set_up()

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

