VariableProcessor \
	--in_data "cmip6_Z_LTDM.nc;cmip6_Z2_LTDM.nc" \
	--out_data "cmip6_Z_mean_stddev.nc" \
	--var "dailymean_zg500,_SQRT(_DIFF(dailymeansq_zg500,_POW(dailymean_zg500,2)))" \
	--varout "dailymean_zg500,stddev_zg500"
