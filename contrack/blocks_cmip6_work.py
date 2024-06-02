####=================================================
# Adapted by Edgar Dolores to  work with ERA5 dataset structure 
# ************nextGEMS project*****************************
# Modified from Daniel Steinfeld 
# For details see:
# Steinfeld, D., 2020: ConTrack - Contour Tracking. GitHub, https://github.com/steidani/ConTrack.
##########


##### FUNCTIONS ######
# ======================================================================
def theta(t,p):
        """
        calculate potential temperature from temperature and pressure
        
        Parameters
        ----------
        p : array_like, 3D or 4D
            Pressure (hPa).
        t : array_like, 3D or 4D
            Temperature (K).
        
        Returns
        -------
        th : array, 3D or 4D
            Theta (potential temperature) (K). Will have same shape
            as p.
        """



        rdcp = 287 / 1003.
        tzero = 273.15
        if np.max(t) < 100:
            t = t + tzero
        th = t*(1000/p)**rdcp
        return(th)

# ======================================================================
# ======================================================================
def portvort(u,v,th,p,lon,lat):
        """
        calculating potential vorticity using fortran subroutine deriv
        
        Need to compile deriv.f90 in console: f2py -c deriv.f90 -m deriv
        
        Parameters
        ----------
        p : array_like, 3D or 4D
            Pressure (hPa).
        th : array_like, 3D or 4D
            potential Temperature (K).
        u,v: array_like, 3D or 4D
            wind (m/s).
            
        
        Returns
        -------
        pv : array, 3D or 4D
            potential vorticity (pvu). Will have same shape as p. 
        """

        #get grid
        nx = len(lons)
        ny = len(lats)
        nz = np.size(p,-3)
        xmin = lon[0]
        ymin = lat[0]
        dx = lons[1] - lons[0]
        dy = np.mean(lats[1:ny] - lats[0:(ny-1)])  # note: small differences between latitudes are ignored
        mdv = -999.99 # missing values

        # CESM data need to be reshaped (TRANSPOSED) for fortran subroutine: (nx, ny, nz) with nz is bottom-up 
        # Steps: squeeze - swap axis - reverse order of z-axis
        u = np.swapaxes(np.squeeze(u), 0, 2)[:,:,::-1]
        v = np.swapaxes(np.squeeze(v), 0, 2)[:,:,::-1]
        th = np.swapaxes(np.squeeze(th), 0, 2)[:,:,::-1]
        p = np.swapaxes(np.squeeze(p), 0, 2)[:,:,::-1]

        # change into SI unit
        p*=100 # hPa to Pa

        #calc PV     
        potvort = deriv.potvort(u, v, th, p, xmin, ymin, dx, dy, mdv)

        #reshape to CESM format
        potvort = np.swapaxes(potvort[:,:,::-1], 0, 2)

       # return data
        return(potvort)

##### ======================================================================
def set_interp2(arg1,p,lon,lat,isobaric): #isobaric ~~ level
        import numpy as np
        import sys; import os


####import fortran subroutine
        sys.path.insert(0, '/home/b/b382006/src/fortran')
        try:
           import extra
        except ImportError:
           msg = """Need to compile deriv with f2py first. By running:
                    f2py -c extra.f90 -m extra
                    Does maybe only work in python 2, or compile it for python 3+.
              """

        ## needed numpy arrays for the interpolation
        dx = 1 
        dy = 1
        xmin = lon[0]
        ymin = lat[0]

        mdv = -999.9

        #### WARNING!!!!!!!!!!!!!!!!
        # sometimes the data need to be reshaped (TRANSPOSED) 
        # for the fortran subroutine:
        # (nx,ny,nz) with nz is bottom-up
        #Steps: squeeze - swap axis - reverse order of z-axis
        pnew  = np.swapaxes(np.squeeze(p),0,2)[:,:,::-1]
        uwind = np.swapaxes(np.squeeze(arg1),0,2)[:,:,::-1]

        iso = np.array(isobaric)
        mz = 1
        #calc
        out1 = extra.int3dm(uwind,pnew,iso,xmin,ymin,dx,dy,mdv)

        #=====reshape to CESM format
        out1 = np.swapaxes(out1,0,2)

        return out1
##### ======================================================================


### FUNCTIONS END ######

### ==========================
### Compute PV AND THETA
#def calc_theta():

#def calc_PV():
### NOT NEEDED for ICON outputs!!!!!!!!!!

### ==================
def calc_VAPV(infile,outdir,outname):
    """
    calculate vertically averaged PV (VAPV) between 500 and 150 hPa with 50hPa steps
    """
    # ======
    #import modules
    import numpy as np
    import netCDF4
    import sys; import os
#    from dypy.small_tools import interpolate

####import fortran subroutine
    sys.path.insert(0, '/home/b/b382006/src/fortran')
    try:
        import extra
    except ImportError:
        msg = """Need to compile deriv with f2py first. By running:
                    f2py -c extra.f90 -m extra
                    Does maybe only work in python 2, or compile it for python 3+.
              """
        print(msg)

   # ======
    #define path
#    pvpath='/work/bm1235/b382006/regridding/icon/from_5km_to_1degree/'

    # name prefix
#    filename = "ngc2009_atm_pl_6h_inst_20200120T000000Z" #ensemble members 001 to 035
    filename = str(infile)
#    outpath = "/work/bm1235/b382006/regridding/icon/from_5km_to_1degree/new_vars"
    outpath = str(outdir)

    fnew = str(outname)
        #=== get PV data
#    infile = pvpath + '/'+ filename +'.nc'
    infile = filename

    con = netCDF4.Dataset(infile, mode='r')

           # load constant variables
#            p0 = con.variables['P0'][:]           # reference pressure [hPa]
#            p0 /= 100   # Pa to hPa
#            hyam = con.variables['hyam'][:]       # Hybrid level coefficient
#            hybm = con.variables['hybm'][:]       # Hybrid level coefficient
    lons = con.variables['longitude'][:]         # longitude
    lats = con.variables['latitude'][:]         # latitude
    levs = con.variables['level'][:]         # pressure
    times = con.variables['time'][:]

        # prepare output array (filled with vapv of each timestep)
    out_vapv = []
   

        # ======
        #loop through timesteps
    n_timesteps = len(times)
    for ii in range(n_timesteps):
                print("{: 4d} of {: 4d}".format(ii + 1, n_timesteps))

            # load time dependent variable
                pv = con.variables['pv'][ii,:]*1000000

### No need of interpolation all levels are given from era5
#            # ======
#            #calculate pressure array: p(k)
#                p = np.zeros_like(pv)
#                for z in range(p.shape[0]):
#                    p[z,...] = levs[z]
#
#            # interpolation and vertically averaging
#                plevel = np.arange(150,550,50)

            #prepare p-grid: not allowed to have somewhere same value as plevels (otherwise interpolate gives nan)
#                for zz in plevel:
#                    index = np.where(p==zz)
#                    for kk in range(len(index[0])):
#                        p[index[0][kk],index[1][kk],index[2][kk]] += 0.001
#
#                pv_int = set_interp2(pv,p,lons,lats,plevel) # ICON from top-to-bottom

                vapv = np.mean(pv,axis=0)

            # add to Out array
                out_vapv.append(vapv)

        #close data
    con.close()

        #from list to array
    out_vapv = np.asarray(out_vapv)

        # ======
        #saving VAPV and time into netcdf file of same format as input netcdf file

    outfile = outpath + '/vapv/' + 'VAPV_' + outname 
    print('saving to file:  \n' + outfile)

    if not os.path.exists(outpath + '/vapv/'):
                os.makedirs(outpath + '/vapv/')

        # create a new NETCDF file
        #open ncfile
    ncfile = netCDF4.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        #create dimensions
    longitude = ncfile.createDimension('longitude', len(lons))
    latitude = ncfile.createDimension('latitude', len(lats))
    time = ncfile.createDimension('time', None)

        #create variables
    longitude = ncfile.createVariable('longitude', np.float32, ('longitude',))
    latitude = ncfile.createVariable('latitude', np.float32, ('latitude',))
    time = ncfile.createVariable('time', np.float32, ('time',))
    pv = ncfile.createVariable('VAPV', np.float32,('time','latitude','longitude'))

       #global attributes
    ncfile.title = "ERA5"
    ncfile.source = "giub@giub.unibe.ch"
    ncfile.comment = "to analyse weather systems "
        #variable attributes
    latitude.units = 'degrees_north'
    longitude.units = 'degrees_east'
    pv.units = 'pvu'
    pv.long_name= 'vertical averaged PV between 500 to 150 hPa'
    time.units = 'hours since 1900-01-01 00:00:00.0'
    time.calendar = "gregorian"
#    time.bounds = "time_bnds"

        #writing data
    longitude[:]=lons
    latitude[:]=lats
    pv[:,:,:] = out_vapv
    time[:] = times
        #close files
    ncfile.close()

###=============================================
def calc_VAPV_clim():
    """
    calculate monthly averaged VAPV 
    Is needed for 
        1) calculating VAPV Anomaly for blocking detection    
    """

    # ======
    #import modules
    from datetime import datetime, timedelta
    import numpy as np
    import netCDF4
    import os

    # ======
    #define path
    inpath = outpath = ''
    present = "b.e112.B20TRLENS.f09_g16.ethz." #ensemble members 001 to 035
    future = "b.e112.BRCP85LENS.f09_g16.ethz." #ensemble members 001 to 035

    #define parameters
    ens_member = [str(item).zfill(3) for item in range(1, 36)]
    years = [str(item) for item in range(1990, 2001)] #present

    # ======
    #set grid dimension and initialize (you can get this information with nz, ny, nx = shape(vapv) )
    nlon = 288
    nlat = 192
    nz = 1
    vapv_clim = []
    ndat = []

    for ii in range(0,12):
        s= np.zeros((nz,nlat,nlon))
        vapv_clim.append(s)
        ndat.append(s)

    # ======
    #loop through years
    for yy in years:
        print(yy)
        # ======
        # get data
        infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + yy + '.nc'
        con = netCDF4.Dataset(infile, mode='r')

        #get variables
        lons = con.variables['longitude'][:]         # longitude
        lats = con.variables['latitude'][:]         # latitude
        times = con.variables['time']
        dtime = [ ]
        dtime = [netCDF4.num2date(xx,units = times.units, calendar=times.calendar) for xx in times]

       #get month
        month = np.asarray([dtime[ii].month for ii in range(0, np.size(dtime))])

        # ======
        #loop through times
        for ii in range(0, np.size(dtime)):
            print(ii)
            # get PV
            pv = con.variables['VAPV'][ii,:]
            npv = con.variables['VAPV'][ii,:]

            npv[npv>-999] = 1; npv[npv<-999] = 0        # 0 where value == -999
            ndat[month[ii]-1] = ndat[month[ii]-1] + npv
            pv[pv<-999] = 0
            vapv_clim[month[ii]-1] = vapv_clim[month[ii]-1] + pv

        #close data
        con.close()

   # ======
    # normalize data
    for ii in range(0,12):
        vapv_clim[ii] = vapv_clim[ii]/(ndat[ii] + 1e-20) #because we cannot devide by zero

    # ======
    #saving monthly VAPV clim into netcdf file of same format as input netcdf file

    outfile = outpath + present + ens_member[0] + '/clim/' + 'VAPV_clim.nc'
    print('saving to file:  \n' + outfile)

    if not os.path.exists(outpath + present + ens_member[0] + '/clim/'):
            os.makedirs(outpath + present + ens_member[0] + '/clim/')

    # create a new NETCDF file
    #open ncfile
    ncfile = netCDF4.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    #create dimensions
    longitude = ncfile.createDimension('longitude', len(lons))
    latitude = ncfile.createDimension('latitude', len(lats))
    time = ncfile.createDimension('time', None)

    #create variables
    longitude = ncfile.createVariable('longitude', np.float32, ('longitude',))
    latitude = ncfile.createVariable('latitude', np.float32, ('latitude',))
    time = ncfile.createVariable('time', np.float32, ('time',))
    pv = ncfile.createVariable('VAPV', np.float32,('time','latitude','longitude'))

    #global attributes
    ncfile.title = "cesm112_LENS"
    ncfile.source = "@met.fu-berlin.de"
    ncfile.comment = "monthly averaged VAPV: to calc blocking with PV index "
    #variable attributes
    latitude.units = 'degrees north'
    longitude.units = 'degrees east'
    pv.units = 'pvu'
    pv.long_name='monthly averaged VAPV from ' + str(years[0]) + ' to ' + str(years[-1])
    time.units = 'months'

    #writing data  
    longitude[:]=lons
    latitude[:]=lats
    pv[:,:,:] = vapv_clim
    time[:] = np.arange(1,13,1)
    #close files
    ncfile.close()

###===================================
def calc_VAPV_anom():
    """
    calculate VAPV anomaly

    """

    # ======
    #import modules
    from datetime import datetime, timedelta
    import numpy as np
    import netCDF4 as nc
    import os

    # ======
    #define path
    inpath = outpath = ''
    present = "b.e112.B20TRLENS.f09_g16.ethz." #ensemble members 001 to 035
    future = "b.e112.BRCP85LENS.f09_g16.ethz." #ensemble members 001 to 035

    #define parameters
    ens_member = [str(item).zfill(3) for item in range(1, 36)]

    # some general settings
    nsmooth = 4  # plus/minus one day running mean for temporal smoothing
    npers = 20   # five day persistence criterion   (CESM: 1 day = 4 timesteps = dt=6h)

    # ======
    # loop over years
    years = np.arange(1990,2001) # present

    for yy in years:
        print(yy)

        # ======
        # get clim data
        climfile = inpath + present + ens_member[0] + '/clim/' + 'VAPV_clim.nc'
        with nc.Dataset(climfile) as ncf:
            clim = ncf.variables['VAPV'][:]
            lons = ncf.variables['longitude'][:]
            lats = ncf.variables['latitude'][:]

       #get instantenous data
        infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy) + '.nc'
        with nc.Dataset(infile) as ncf:
            vapv = ncf.variables['VAPV'][:]
            times = ncf.variables['time'] #  days since refdate
            date = [nc.num2date(xx,units = times.units, calendar=times.calendar) for xx in times]
            month = np.asarray([date[ii].month for ii in range(0, np.size(date))])

        #taking data availability into account in previous/next year (need 5 days in each)
        if (yy != years[-1]):
            infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy+1) + '.nc'
            with nc.Dataset(infile) as ncf:
                vapv_past = ncf.variables['VAPV'][:npers]
                vapv = np.concatenate((vapv,vapv_past), axis = 0)
                month = np.concatenate((month,np.ones(npers)), axis = 0)

        if (yy != years[0]):
            infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy-1) + '.nc'
            with nc.Dataset(infile) as ncf:
                vapv_pre = ncf.variables['VAPV'][-npers:]
                vapv = np.concatenate((vapv_pre,vapv), axis = 0)
                month = np.concatenate((np.full(npers, 12),month), axis = 0)

        # ======
        # Note that index start at 0: month-1 to get right month in clim
        del times,date
        month -= 1
        month = month.astype(int)
        ndat = len(month)
        # calc anomaly from monthly avarage
        anom = [vapv[x,...] - clim[month[x],...] for x in range(0,len(month))]
        anom = np.asarray(anom)
        # SH reverse sign
        anom[:,:96,:] = anom[:,:96,:] * (-1)

        # ======
        # do temporal smoothing (2 day running mean: 2 *nsmooth)
        anom_save = anom[0:(0+nsmooth-1),...]
        for ii in range(0,ndat-2*nsmooth):
            for jj in range(0,2*nsmooth):
                anom[ii,...] = anom[ii,...] + anom[ii+jj,...]
            anom[ii,...] = anom[ii,...] / (2*nsmooth+1)
        anom[(0+nsmooth):(ndat-nsmooth),...] = anom[0:(ndat-2*nsmooth),...]
        anom[0:(0+nsmooth-1),...] = anom_save

        # set lat[0] and lat[-1] (North and South Pole) to zero, so that block cannot jump between Hemis (happens sometimes...)
        anom[:,0,:] = 0
        anom[:,-1,:] = 0

        # ======
        #get time
        infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy) + '.nc'
        with nc.Dataset(infile) as ncf:
            times = ncf.variables['time'][:] #  days since refdate
        if (yy != 2000):
            infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy+1) + '.nc'
            with nc.Dataset(infile) as ncf:
                times_past = ncf.variables['time'][:20]
                times = np.concatenate((times, times_past), axis = 0)
        if (yy != 1990):
            infile = inpath + present + ens_member[0] + '/vapv/' + 'VAPV_' + str(yy-1) + '.nc'
            with nc.Dataset(infile) as ncf:
                times_pre = ncf.variables['time'][-20:]
                times = np.concatenate((times_pre, times), axis = 0)

        # ======
       #saving VAPV anom into netcdf file of same format as input netcdf file

        outfile = outpath + present + ens_member[0] + '/anom/' + 'Vanom' + str(yy)
        print('saving to file:  \n' + outfile)

        if not os.path.exists(outpath + present + ens_member[0] + '/anom/'):
                os.makedirs(outpath + present + ens_member[0] + '/anom/')

        # create a new NETCDF file    
        #open ncfile
        ncfile = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')

        #create dimensions
        longitude = ncfile.createDimension('longitude', len(lons))
        latitude = ncfile.createDimension('latitude', len(lats))
        time = ncfile.createDimension('time', None)

        #create variables
        longitude = ncfile.createVariable('longitude', np.float32, ('longitude',))
        latitude = ncfile.createVariable('latitude', np.float32, ('latitude',))
        time = ncfile.createVariable('time', np.float32, ('time',))
        VAPVanom = ncfile.createVariable('VAPVanom', np.float32,('time','latitude','longitude'))

        #global attributes
        ncfile.title = "blocks in cesm112_LENS"
        ncfile.source = "@met.fu-berlin.de"
        ncfile.comment = "to analyse weather systems "
        #variable attributes
        latitude.units = 'degrees_north'
        longitude.units = 'degrees_east'
        VAPVanom.long_name= 'vertical averaged PV anomaly between 500 to 150 hPa'
        VAPVanom.units = 'PVU'
        #pv.missing_value = -999
        time.units = 'days since 1850-01-01 00:00:00'
        time.calendar = "noleap"
        time.bounds = "time_bnds"

        #writing data  
        time[:] = times
        latitude[:]=lats
        longitude[:]=lons
        VAPVanom[:,:,:] = anom

        #close files
        ncfile.close()

####===========================
####***********************************
####===========================

def blocks_freq():
    """ calculate blocking frequency for the CESM: annual and for each season
        can be visualized with plot_block_freq
    """
    # =====
    # load modules
    import numpy as np
    import netCDF4 as nc
    from datetime import datetime, timedelta
    import pickle

    # =====
    # define path
    inpath = outpath = ''
    present = "b.e112.B20TRLENS.f09_g16.ethz." #ensemble members 001 to 035
    future = "b.e112.BRCP85LENS.f09_g16.ethz." #ensemble members 001 to 035

    #define parameters
    ens_member = [str(item).zfill(3) for item in range(1, 36)]

    # define grid
    dx = 1.25   # taken from inputfile
    dy = 0.94240837696335078 # taken from inputfile
    lon = np.arange(0,360,dx)
    lat = np.arange(-90,90.5,dy)
    nlon = len(lon)
    nlat = len(lat)
    #clean up
    del dx, dy, lon, lat


   # =====
    # initialize array for annual and seasonal frequencies
    block_sf = []
    block_f = np.zeros((nlat,nlon))
    for ii in range(0,4): # for season use 0:4
        s= np.zeros((nlat,nlon))
        block_sf.append(s)
    nd_seas = np.zeros((4,1))
    nd = 0

    # =====
    # loop over year
    for yy in range(1990,2001):
        print(yy)

        # =====
        #get data
        # PV index
        infile = inpath + present + ens_member[0] + '/block/' + 'BLOCKS' + str(yy) + '.nc'
        # AGP index
        #infile = inpath + present + ens_member[0] + '/agp/' + 'BLOCKSagp' + str(yy) + '.nc'

        with nc.Dataset(infile) as ncf:
            #print(ncf)
            flag = ncf.variables['FLAG'][:]
            # convert hour to datetime
            times = ncf.variables['time'] #  days since refdate
            date = [nc.num2date(xx,units = 'days since 1850-01-01 00:00:00', calendar="noleap") for xx in times]

       # =====
        # need only zeros and ones
        flag[flag>1] = 1

        # get dates for each season
        month = np.asarray([date[ii].month for ii in range(0, np.size(date))])
        year = np.asarray([date[ii].year for ii in range(0, np.size(date))])
        ind_seas = []
        ind_seas.append(np.where((month<=2) | (month==12) & (year==yy))) #DJF
        ind_seas.append(np.where((month>=3) & (month<=5) & (year==yy))) #MAM
        ind_seas.append(np.where((month>=6) & (month<=8) & (year==yy))) #JJA
        ind_seas.append(np.where((month>=9) & (month<=11) & (year==yy))) # SON
        #clean up
        del month, year, times

       #sum over season
        for jj in range(0,4):
            nd_seas[jj] = nd_seas[jj] + len(ind_seas[jj][0])
            for ii in ind_seas[jj][0]:
                block_sf[jj] = block_sf[jj] + flag[ii,...]
        #sum over whole year
        for ii in range(0,len(date)):
            block_f = block_f + flag[ii,...]
        nd = nd + len(date)

        #clean up
        del flag, date, ind_seas

    # normalize
    for jj in range(0,4):
        block_sf[jj] = block_sf[jj] / nd_seas[jj]
    block_f = block_f / nd

    # Saving the objects in file:
    with open(outpath + present + ens_member[0] + '/data/' + 'blocks_freq_seas', 'wb') as f:
        pickle.dump([block_f, block_sf], f)



def calc_gh_anom():
    """
    calculate gh anomaly

    """

    # ======
    #import modules

    from contrack import contrack
    import xarray as xr
    import datetime
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import numpy as np
    from matplotlib import cm
    import sys, os, argparse
    import bottleneck


    ### open cmip6 files
    xr_in=xr.open_mfdataset('/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/factory/cmip6_zg500_*')
    var="zg500"
    prefix="MPI-ESM1-2-LR"
    outdir="/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/factory/z500/"
    print('start preprocessing ...')

    ### compute climatology
    window=31
    groupby="dayofyear"
    clim = xr_in[var].groupby("time.dayofyear").mean("time")
    clim2 = clim.rolling(**{groupby:window}, center=True).mean().fillna(clim[-window:].mean(dim=groupby))

    ### compute anomaly 
    smooth=8
    groupby='dayofyear'
    anom_gh=(xr_in[var].groupby('time.dayofyear') - clim2).rolling(time=smooth, center=True).mean()

    ### cts for computing geopotential height
#    g = 9.80665  # m s**-2
#    anom_gh=anom/g
#    ### attributes
#    anom_gh.attrs['units'] = 'm'
#    anom_gh.attrs['long_name']= 'Geopotential Height Anomaly'


    ### saving by year and adding 5 days in december
    y0=anom_gh.time[0].dt.year
    yn=anom_gh.time[-1].dt.year
    npers=20 ###5 days
    for yy in range(int(y0),int(yn)+1):
        anoms = anom_gh.sel(time=anom_gh.time.dt.year.isin(yy))
 #      if (yy != yn):
 #       anom_tmp = anom_gh.sel(time=anom_gh.time.dt.year.isin(yy+1))
 #       anom_past = anom_tmp[:npers]
 #       anoms = xr.concat([anoms, anom_past], dim="time")
 #
 #      if (yy != y0):
 #       anom_tmp = anom_gh.sel(time=anom_gh.time.dt.year.isin(yy-1))
 #       anom_pre = anom_tmp[-npers:]
 #       anoms = xr.concat([anom_pre,anoms], dim = "time")

    #print(anoms)
        outfile=prefix+"-anom-z500_"+str(yy)+".nc"
        print(outfile)

        anoms.to_netcdf(outdir+outfile)



# ======================================================================================================================================

###### SELECT FUNTION TO RUN ###############
print("calc_VAPV()")

import sys
#infile = sys.argv[1]
#outdir = sys.argv[2]
#outname = sys.argv[3]
#calc_VAPV(infile,outdir,outname)

print("calc_gh_anom()")
calc_gh_anom()

############################################
# ==========================================================================================================================================

