### blocking method based on the original Code of Daniel Steinfeld 
from contrack import contrack
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import sys

inpath = sys.argv[1]
outdir = sys.argv[2]
outname = sys.argv[3]

print(inpath)

####====INPUT
###===ngc2009
#xr_in=xr.open_mfdataset('/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree/pl/ngc2009_atm_pl_6h_inst_2*')
###===ngc2013
#xr_in=xr.open_mfdataset(inpath+'/ngc2013_atm_pl_3h_inst_204*',concat_dim='time',combine='nested')
###===ngc2013 time missed solve???
xr_in=xr.open_dataset('/scratch/b/b382006/nextgems/cycle2/ngc2013/pl/ngc2013_z500_6h.nc')

###===OUTPUT
###===ngc2009
#outpath = '/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree/block/ANO_Z500/'
###==ngc2013
#outpath = '/work/bm1235/b382006/regridding/icon/ngc2013_10km_to_1degree/block/ANO_Z500/'
outpath = '/scratch/b/b382006/nextgems/cycle2/ngc2013/block/ANO_Z500/'

###==subfix ngc2009
#subfix='-icon-ngc2009'

###==subfix ngc2013
subfix='-icon-ngc2013-p1'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'
outfile_table='BLOCKS'+subfix+'.csv'

###NO APPLY anymore this:
###resample *****ngc2013******* frequency = 3 hourly
###xr_in=xr_in.resample(time='6H').mean(dim='time',keep_attrs=True) #6 hourly mean

## rename zfull
xr_in=xr_in.rename({'zfull':'z_height'})

###select some years ###== new folder
years=[2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030]
xr_in=xr_in.z_height.sel(time=xr_in.z_height.time.dt.year.isin([years]))

### Take geopotential
#in_z_500=xr_in.z_height.sel(plev=50000.0)
in_z_500=xr_in.sel(plev=50000.0) ### new folder
z_dataset=xr.DataArray.to_dataset(in_z_500)


### intitate
block = contrack()
block.read_xarray(z_dataset)

# calculate geopotential height
#block.calculate_gph_from_gp(gp_name='z',
#                            gp_unit='m**2 s**-2',
#                            gph_name='z_height')

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
###block


# Out[] Xarray dataset with 2707 time steps.
#            Available fields: z, z_height, anom, flag

# Hint: In case you want to use a more objective threshold, e.g., the 90th percentile of the Z500 anomaly winter distribution over 50°-80°N, do:
# threshold = block['anom'].sel(latitude=slice(80, 50)).quantile([0.90], dim='time').mean() # 177gmp

# save to disk
block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

###==== flag = output of block.run_contrack(), variable = input variable to calculate intensity and center of mass
block_df = block.run_lifecycle(flag='flag', variable='anom')
block_df.to_csv(outpath+'/'+outfile_table, index=False)
