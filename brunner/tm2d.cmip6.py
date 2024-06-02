import sys
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import cm
import os

sys.path.append("/home/b/b382006/tools/Blocking/brunner/blocking/")
from BlockingDetection import Blocking

### input
inpath="/scratch/b/b382006/cmip6/"
model="MPI-ESM1-2-LR/"
member_id="r9i1p1f1" 
#r10i1p1f1, r2i1p1f1,  r4i1p1f1,  r6i1p1f1,  r8i1p1f1,  
#r1i1p1f1,  r3i1p1f1,  r5i1p1f1,  r7i1p1f1,  r9i1p1f1, 
y0=1979
yn=2015
years=np.arange(y0,yn,1)

for year in years:
    print(year)
    ### Compute yearly
#    infile=inpath+model+member_id+"/daily/cmip6_z500_day_"+str(year)+".nc"
    infile=inpath+model+member_id+"/daily/cmip6_zg500_"+str(year)+"_day.nc"

###==subfix ngc2013
#subfix='-ifs-daily'

###==output
    dir_res=inpath+model+member_id+"/block/ABS_500"
    ofile= "block_ABS-Z500_"+str(year)+".nc"

    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    OUTPATH=dir_res+'/'+ofile

    print(infile,OUTPATH)

    data=xr.open_dataset(infile)


###=== rename zg
#    data=data.rename({'z':'GeopotentialHeight'})
    data=data.rename({'zg500':'GeopotentialHeight'})


###=== intitate
    blk = Blocking()

    blk.import_xarray(data)

#blk.calculate_gph_from_gp() # calculate geopotential height

    blk.set_up(time_name='time',longitude_name='lon', latitude_name='lat')

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

    data.close()
