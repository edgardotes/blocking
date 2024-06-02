### blocking method based on the original Code of Daniel Steinfeld 
from contrack import contrack
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys


for yy in range(2020,2025):

####====INPUT
###===ifs 4.4 fesom 5
    print("1. Load files")
    
###===IFS_4.4-FESOM_5
#    xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/z500/ifs_c3-anom-z500_'+str(yy)+'.nc')
###===IFS_28-NEMO_25
    xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_28-NEMO_25/z500/ifs-c3-anom-z500_'+str(yy)+'.nc')

###===OUTPUT
###==ifs 4.4 km
#    outpath = '/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/block/ANO_Z500/'
###===IFS_28-NEMO_25
    outpath = '/scratch/b/b382006/nextgems/cycle3/IFS_28-NEMO_25/block/ANO_Z500/'

###==subfix ifs 4.4
#    subfix='-ifs-c3-'+str(yy)+''
###==subfix ifs 28 
    subfix='-ifs-c3-'+str(yy)+''

###OUTFILE NAME
    outfile_flag='BLOCKS'+subfix+'.nc'
    outfile_table='BLOCKS'+subfix+'.csv'

###varname
    var="z" 

    print('start preprocessing ...')

### intitate
    block = contrack()
    block.read_xarray(xr_in)

# Hint: Use block.set_up(...) to do consistency check and set (automatically or manually) names of dimension ('time', 'latitude', 'longitude')
    block.set_up(force=True)
    block.ds=block.ds.compute()


# Finally, track blocking anticyclones (>=150gmp, 50% overlap twosided, 5 timesteps persistence (here 5 days))
    block.run_contrack(variable=var,
                   threshold=160,
                   gorl='>=',
                   overlap=0.7,
                   persistence=20,
                   twosided=True)

# output: variable 'flag'. 440 blocking systems tracked. Each blocking system is identified by a unique flag/ID.
###block


# Out[] Xarray dataset with 2707 time steps.
#            Available fields: z, z_height, anom, flag

# Hint: In case you want to use a more objective threshold, e.g., the 90th percentile of the Z500 anomaly winter distribution over 50°-80°N, do:
# threshold = block['anom'].sel(latitude=slice(80, 50)).quantile([0.90], dim='time').mean() # 177gmp

# save to disk
    block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

###==== flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
    block_df = block.run_lifecycle(flag='flag', variable=var)
    block_df.to_csv(outpath+'/'+outfile_table, index=False)
