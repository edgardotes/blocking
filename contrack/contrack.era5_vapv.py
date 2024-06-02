### blocking method based on the original Code of Daniel Steinfeld 
from contrack import contrack
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm


#### FUNCTIONS #####
def print_now():
    print('time:'+str(datetime.datetime.now()))


####====INPUT
xr_in=xr.open_mfdataset('/work/bm1235/b382006/era5/vapv/VAPV_era5_regio_prl_202*',concat_dim='time',combine='nested')

###===OUTPUT
outpath = '/work/bm1235/b382006/era5/block/ANO_VAPV'

###==subfix era5
subfix='-era5-ano_2020-2022'


###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'
outfile_table='BLOCKS'+subfix+'.csv'


print('start preprocessing ...')

###change to yyyy mm dd
xr_in['time']=xr_in.indexes['time'].normalize()

###Drop dimension time_bnds
#dxr_in=xr_in.drop_dims('bnds')


### Take geopotential
#in_z_500=xr_in.z.sel(plev=50000.0)
#z_dataset=xr.DataArray.to_dataset(in_z_500)

#z_dataset=xr_in ##taking directly Z


### intitate
block = contrack()
block.read_xarray(xr_in)


# Hint: Use block.set_up(...) to do consistency check and set (automatically or manually) names of dimension ('time', 'latitude', 'longitude')
block.set_up(force=True)
block.ds=block.ds.compute()

### compute anomalies
print_now()
print('computing vapv clim/anom')
block.calc_anom(variable='VAPV',
        smooth=8,
        window=31,
        groupby='dayofyear')

####Identify and track blocks
block.run_contrack(variable='anom',
                   threshold=-1,
                   gorl='<=',
                   overlap=0.7,
                   persistence=20,
                   twosided=True)

# save to disk
block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

###==== flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
block_df = block.run_lifecycle(flag='flag', variable='anom')
block_df.to_csv(outpath+'/'+outfile_table, index=False)
