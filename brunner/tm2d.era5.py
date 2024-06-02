###=== blocking method based on the original Code of Lukas Brunner
#### WORKING WITH ERA5
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys
sys.path.append("/home/b/b382006/tools/Blocking/brunner/blocking/")
from BlockingDetection import Blocking

inpath = sys.argv[1]
outdir = sys.argv[2]
outname = sys.argv[3]

print(inpath,outname)

####====INPUT
##xr_in=xr.open_mfdataset(inpath+'/era5_regio_prl_*',concat_dim='time',combine='nested')
xr_in=xr.open_dataset(inpath+'/'+outname,chunks={'time': 1})

print('data loaded ...')

###===OUTPUT
outpath = '/work/bm1235/b382006/era5/block/ABS_Z500/'

###==subfix ngc2013
subfix='-era5_1959-1989'

###OUTFILE NAME
outfile_flag1='IB_BLOCKS'+subfix+'.nc'
outfile_flag2='LSB_BLOCKS'+subfix+'.nc'
outfile_flag3='BLOCKS'+subfix+'.nc'


print('start preprocessing ...')

###change to yyyy mm dd
xr_in['time']=xr_in.indexes['time'].normalize()

###Drop dimension time_bnds
dxr_in=xr_in.drop_dims('bnds')

### rename zfull
#xr_in=xr_in.rename({'zfull':'GeopotentialHeight'})

### Take geopotential at 500 hPa
#in_z_500=xr_in.GeopotentialHeight.sel(plev=50000.0)

###daily mean
#dataset_z=xr_in.resample(time='1D').mean(dim='time',keep_attrs=True) 
dataset_z=dxr_in


###select months
###pd_z = daily_z.sel(time=daily_z.time.dt.month.isin([1, 2, 12]))
###dataset_z=xr.DataArray.to_dataset(pd_z)
#dataset_z=xr.DataArray.to_dataset(daily_z)

### configure attributes 
dataset_z.time.attrs['standard_name'] = 'time'
dataset_z.time.attrs['axis'] = 'T'
#print(dataset_z)


### Outpath line
###OUTPATH = dataset_z.replace(
###    'geopotential', 'blocking').replace('.nc', '_bf_3D.nc')
OUTPATH1=outpath+'/'+outfile_flag1
OUTPATH2=outpath+'/'+outfile_flag2
OUTPATH3=outpath+'/'+outfile_flag3

print('start main program ...')
###=== intitate
blk = Blocking()
blk.import_xarray(dataset_z)

###=== pre-processing
###blk.get_time_subset(months='DJF')
###blk.calculate_daily_mean()

blk.calculate_gph_from_gp() # calculate geopotential height
#blk.ds.time

####configurate and run
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

print('Saving blockings ...')

###=== save to disk 
#blk.save(OUTPATH1, 'IB')
#blk.save(OUTPATH2, 'ExtendedIB')
blk.save(OUTPATH3, 'Blocking')
