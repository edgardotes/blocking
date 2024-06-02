DetectBlobs \
	--in_data_list cmip6_DB_files.txt \
	--out_list cmip6_blocktag_files.txt \
	--thresholdcmd "_DIFF(zg500,threshold_zg500),>=,0,0" \
	--minabslat 25 \
	--maxabslat 75 \
	--geofiltercmd "area,>,1e6km2" \
	--tagvar "block_tag"
