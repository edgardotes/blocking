FourierFilter \
	--in_data cmip6_Z_mean_stddev.nc \
	--out_data cmip6_Z_mean_stddev_timesmoothed.nc \
	--var "dailymean_zg500,stddev_zg500" \
	--dim "time" \
	--modes 4

FourierFilter \
	--in_data cmip6_Z_mean_stddev_timesmoothed.nc \
	--out_data cmip6_Z_mean_stddev_smoothed.nc \
	--var "stddev_zg500" \
	--preserve "dailymean_zg500" \
	--dim "lon" --modes 2
