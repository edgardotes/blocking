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

### blocking method based on the original Code of Daniel Steinfeld 
xr_in=xr.open_mfdataset('/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree/vapv/VAPV_ngc2009*',concat_dim='time',combine='nested')


outpath = '/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree/block/ANO_VAPV/'

subfix='-icon-ngc2009'
outfile_flag='BLOCKS'+subfix+'.nc'
outfile_table='BLOCKS'+subfix+'.csv'

### intitate
block = contrack()
block.read_xarray(xr_in)


block.set_up(force=True)
block.ds = block.ds.compute() ### look at this!

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

# flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
block_df = block.run_lifecycle(flag='flag', variable='anom')
block_df.to_csv(outpath+'/'+outfile_table, index=False)
