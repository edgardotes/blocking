#%pip install contrack
###== load contrack
from contrack import contrack
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm

#ifs_data=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/z500/z500_ifs_c3.nc')
ifs_data=xr.open_dataset('/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/z500/z500_ifs_c3_test_v3.nc')

outpath = '/scratch/b/b382006/nextgems/cycle3/IFS_4.4-FESOM_5/block/ANO_Z500'

#subfix='-ifs_4.4-fesom_5'
subfix='-ifs_4.4-fesom_5_full_160-07'

###OUTFILE NAME
outfile_flag='BLOCKS'+subfix+'.nc'
outfile_table='BLOCKS'+subfix+'.csv'

### intitate
block = contrack()
#block.read_xarray(tmp)
block.read_xarray(ifs_data)

# calculate geopotential height
block.calculate_gph_from_gp(gp_name='z',
                            gp_unit='m**2 s**-2',
                            gph_name='z_height')

block.set_up(force=True)
block.ds=block.ds.compute()

# calculate Z500 anomaly (temporally smoothed with a 2 d running mean) with respect to the 31-day running mean (long-term: 30 years) climatology
block.calc_anom(variable='z_height',
                smooth=8,
                window=31,
                groupby='dayofyear')

# Finally, track blocking anticyclones (>=150gmp, 50% overlap twosided, 5 timesteps persistence (here 5 days))
block.run_contrack(variable='anom',
                   threshold = 160,
                   gorl='>=',
                   overlap=0.7,
                   persistence=20,
                   twosided=True)

# Hint: In case you want to use a more objective threshold, e.g., the 90th percentile of the Z500 anomaly winter distribution over 50°-80°N, do:
#threshold = block['anom'].sel(latitude=slice(80, 50)).quantile([0.90], dim='time').mean() # 177gmp


# save to disk
block['flag'].to_netcdf(outpath+'/'+outfile_flag,unlimited_dims='time')

block_df = block.run_lifecycle(flag='flag', variable='anom')
block_df.to_csv(outpath+'/'+outfile_table, index=False)



