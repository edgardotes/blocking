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

### One of these models
#models=["HadGEM3-GC31-MM","CNRM-ESM2-1","UKESM1-0-LL",
#        "HadGEM3-GC31-LL","MRI-ESM2-0",
#        "CNRM-CM6-1","ACCESS-CM2",
#        "EC-Earth3","BCC-CSM2-MR",
#        "MPI-ESM1-2-HR", "CESM2-WACCM",
#        "MIROC6","MPI-ESM1-2-LR","CESM2"]

### input
inpath="/scratch/b/b382006/cmip6/"
model="ACCESS-CM2"
member_id="r1i1p1f1" 
#r10i1p1f1, r2i1p1f1,  r4i1p1f1,  r6i1p1f1,  r8i1p1f1,  
#r1i1p1f1,  r3i1p1f1,  r5i1p1f1,  r7i1p1f1,  r9i1p1f1,
#r11i1p1f1 
y0=1979 #1979
yn=2014 #2015
years=np.arange(y0,yn,1)

###full dataset
infile=inpath+model+"/"+member_id+"/daily/post_cmip6_Z50000_day.nc"
print(infile)
xr_data=xr.open_dataset(infile)
xr_data=xr_data.isel(plev=0,drop=True) # remove plev

#xr_data['time'].attrs['calendar']='noleap'
#print(xr_data['time'])

#print(xr_data)

for year in years:
###==output
    dir_res=inpath+model+"/"+member_id+"/block/ABS_500"
    ofile= "block_ABS-Z500_"+str(year)+".nc"

    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    OUTPATH=dir_res+'/'+ofile

    print(infile,OUTPATH)

    ### Compute yearly
    data=xr_data.sel(time=xr_data.time.dt.year.isin([year]))

###=== intitate
    blk = Blocking()

    blk.import_xarray(data)

    blk.calculate_gph_from_gp(gp_name='Z') # calculate geopotential height

#    print(data.time[0:365])
#    print(data.time[364])
#    delta = np.unique((data.time[1:] - data.time[:-1]).days)
#    print('delta',delta)
#    var=data.time.to_index()
#    print(var)
#    delta = np.unique((var[1:] - var[:-1]).astype('timedelta64[D]'))
#    delta = (var[1:] - var[:-1]).astype('timedelta64[D]')
#    print('delta',delta[58])


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
    blk
    blk.save(OUTPATH, 'Blocking')

    data.close()
