### blocking method based on the original Code of Daniel Steinfeld 
from contrack import contrack
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm

y0=1979
yn=2015
var="zg500"

####====INPUT
for yy in range(int(y0),int(yn)+1):
    print("1. Load files")
    xr_in=xr.open_dataset("/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/factory/z500/MPI-ESM1-2-LR-anom-z500_"+str(yy)+".nc")

###===OUTPUT
    outpath="/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/block/ANO_Z500"

###==subfix era5
    subfix='-MPI-ESM1-2-LR-ano_6h-'+str(yy)+''


###OUTFILE NAME
    outfile_flag='BLOCKS'+subfix+'.nc'
    outfile_table='BLOCKS'+subfix+'.csv'


    print('start preprocessing ...')


### intitate
    block = contrack()
    block.read_xarray(xr_in)

# Hint: Use block.set_up(...) to do consistency check and set (automatically or manually) names of dimension ('time', 'latitude', 'longitude')

    block.set_up(force=True)
    block.ds=block.ds.compute()

# Hint: Use 'clim=...' to point towards an existing climatological mean (useful for weather forecasts)
# output: variable 'anom'.

# Finally, track blocking anticyclones (>=160gmp, 70% overlap twosided, 20 timesteps persistence (here 5 days))
    block.run_contrack(variable=var,
                   threshold=160,
                   gorl='>=',
                   overlap=0.7,
                   persistence=20,
                   twosided=True)

# output: variable 'flag'. 440 blocking systems tracked. Each blocking system is identified by a unique flag/ID.

## Hint: In case you want to use a more objective threshold, e.g., the 90th percentile of the Z500 anomaly winter distribution over 50°-80°N, do:
## threshold = block['anom'].sel(latitude=slice(80, 50)).quantile([0.90], dim='time').mean() # 177gmp

# save to disk
    block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

###==== flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
    block_df = block.run_lifecycle(flag='flag', variable=var)
    block_df.to_csv(outpath+'/'+outfile_table, index=False)

    xr_in.close()
