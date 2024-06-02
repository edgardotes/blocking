StitchBlobs \
	--in_list cmip6_blocktag_files.txt \
	--out_list cmip6_blockid_files.txt \
	--var "block_tag" \
	--mintime "5d" \
	--min_overlap_prev 20 \
	--flatten
