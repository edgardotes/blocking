VariableProcessor \
	--in_data cmip6_Z_mean_stddev_smoothed.nc \
	--out_data cmip6_threshold_Z_filtered.nc \
	--var "_SUM(dailymean_zg500,_MAX(100.0,_PROD(1.5,stddev_zg500)))" \
	--varout "threshold_zg500"
