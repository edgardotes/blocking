from contrack import contrack
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm


### blocking method based on the original Code of Daniel Steinfeld 
xr_in=xr.open_mfdataset('/work/bm1235/b382006/regridding/ifs/tco1279-orca025_9km_to_1degree/Y*/ICMUAhr2n+0*')

###===OUTPUT
outpath = '/work/bm1235/b382006/regridding/ifs/tco1279-orca025_9km_to_1degree/block/ANO_Z500/'

###==subfix ngc2013
subfix='-ifs-tco1279'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'
outfile_table='BLOCKS'+subfix+'.csv'

### Take geopotential
in_z_500=xr_in.z.sel(plev=50000.0)
z_dataset=xr.DataArray.to_dataset(in_z_500)


### intitate
block = contrack()
block.read_xarray(z_dataset)

# calculate geopotential height
block.calculate_gph_from_gp(gp_name='z',
                            gp_unit='m**2 s**-2',
                            gph_name='z_height')

# Hint: Use block.set_up(...) to do consistency check and set (automatically or manually) names of dimension ('time', 'latitude', 'longitude')

block.set_up(force=True)
block.ds=block.ds.compute()

# calculate Z500 anomaly (temporally smoothed with a 2 d running mean) with respect to the 31-day running mean (long-term: 30 years) climatology
block.calc_anom(variable='z_height',
                smooth=8,
                window=31,
                groupby='dayofyear')

# Hint: Use 'clim=...' to point towards an existing climatological mean (useful for weather forecasts)
# output: variable 'anom'.

# Finally, track blocking anticyclones (>=150gmp, 50% overlap twosided, 5 timesteps persistence (here 5 days))
block.run_contrack(variable='anom',
                   threshold=160,
                   gorl='>=',
                   overlap=0.7,
                   persistence=20,
                   twosided=True)

# output: variable 'flag'. 440 blocking systems tracked. Each blocking system is identified by a unique flag/ID.
#block


# Hint: In case you want to use a more objective threshold, e.g., the 90th percentile of the Z500 anomaly winter distribution over 50°-80°N, do:
# threshold = block['anom'].sel(latitude=slice(80, 50)).quantile([0.90], dim='time').mean() # 177gmp

# save to disk
block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

###==== flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
block_df = block.run_lifecycle(flag='flag', variable='anom')
block_df.to_csv(outpath+'/'+outfile_table, index=False)
