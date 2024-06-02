#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculates a blocking index following Brunner et al. (2017)
Data 0.5 degree
"""
import sys
sys.path.append("/Users/edgardolorestesillos/Documents/NEXTGEMS/tools/blocking-master/blocking/")
from BlockingDetection import Blocking

INPATH = '/Users/edgardolorestesillos/Documents/NEXTGEMS/tools/blocking-master/data/geop_1degree.nc'
#INPATH = '/Users/edgardolorestesillos/Documents/NEXTGEMS/tools/blocking-master/data/geop_05degree.nc'
OUTPATH = INPATH.replace(
    'geopotential', 'blocking').replace('.nc', '_bf_3D.nc')

blk = Blocking()
blk.read(INPATH)

#blk.get_time_subset()
#blk.calculate_daily_mean()
blk.calculate_gph_from_gp()

blk.set_up()

blk.calculate_gradients(delta_degree=15)


blk.calculate_ib(
    gradient_equator_below=0,
    gradient_pole_below=-10,
    gradient_equator2_above=5)


blk.calculate_eib(min_extent_degree=15)

#longitude_pm_degree, latitude_pm_degree (float, optional):
#          Plus/minus degree longitude and latitude around the center that
#          the block is allowed to move during the given time range
#          (default=7.5/2.5 lon/lat)

blk.calculate_blocking(
    stationary_pm_days=2,
    longitude_pm_degree=8,
    latitude_pm_degree=2)

#blk.reduce_to_1D([0, 75])

blk.save(OUTPATH, 'Blocking')
