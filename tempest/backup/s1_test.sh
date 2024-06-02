Climatology \
	--in_data_list cmip6_Z_list.txt \
	--out_data cmip6_Z_LTDM.nc \
	--var "zg500" \
	--period "daily" \
	--type "mean" \
	--missingdata

Climatology \
	--in_data_list cmip6_Z_list.txt \
	--out_data cmip6_Z2_LTDM.nc \
	--var "zg500" \
	--period "daily" \
	--type "meansq" \
	--missingdata
