#!/bin/bash
#SBATCH --job-name=run_block_rev
#SBATCH -p shared
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -o run_block.out
#SBATCH -e run_block.err
#SBATCH --mem=20Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
module load python3
module load cdo netcdf-c nco emacs

#==extract z500 and compute mean day
loadpath=/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/cmip6_zg500_*
savepath=/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/daily
factory=/scratch/b/b382006/cmip6/MPI-ESM1-2-LR/r2i1p1f1/factory/tempest

#### create folder to save data
if [ ! -d ${savepath}  ]; then
   mkdir -p ${savepath};
fi

#cdo -b F32 daymean -mergetime $loadpath $savepath/cmip6_z500_day.nc

###=== prepare for main program
#python rename_dims_cmip6.py

###===start running tempest routines
source activate tempest_env

###==== Compute z500 climatology

###create lists
#ls ${savepath}/post_* > cmip6_z500_list.txt
#
#Climatology \
#        --in_data_list cmip6_z500_list.txt \
#        --out_data ${factory}/cmip6_z500_LTDM.nc \
#        --var "z" \
#        --period "daily" \
#        --type "mean" \
#        --missingdata
#
#Climatology \
#        --in_data_list cmip6_z500_list.txt \
#        --out_data ${factory}/cmip6_z500-2_LTDM.nc \
#        --var "z" \
#        --period "daily" \
#        --type "meansq" \
#        --missingdata
#
###== calculating standard deviation
#
#VariableProcessor \
#        --in_data "${factory}/cmip6_z500_LTDM.nc;${factory}/cmip6_z500-2_LTDM.nc" \
#        --out_data "${factory}/cmip6_z500_mean_stddev.nc" \
#        --var "dailymean_z,_SQRT(_DIFF(dailymeansq_z,_POW(dailymean_z,2)))" \
#        --varout "dailymean_z,stddev_z"
#
####=== 4-mode  fourier filter
#FourierFilter \
#        --in_data ${factory}/cmip6_z500_mean_stddev.nc \
#        --out_data ${factory}/cmip6_z500_mean_stddev_timesmoothed.nc \
#        --var "dailymean_z,stddev_z" \
#        --dim "time" \
#        --modes 4
#
#FourierFilter \
#        --in_data ${factory}/cmip6_z500_mean_stddev_timesmoothed.nc \
#        --out_data ${factory}/cmip6_z500_mean_stddev_smoothed.nc \
#        --var "stddev_z" \
#        --preserve "dailymean_z" \
#       --dim "lon" --modes 2
#
###=== computing thereshold
#VariableProcessor \
#        --in_data ${factory}/cmip6_z500_mean_stddev_smoothed.nc \
#        --out_data ${factory}/cmip6_threshold_z500_filtered.nc \
#        --var "_SUM(dailymean_z,_MAX(100.0,_PROD(1.5,stddev_z)))" \
#        --varout "threshold_z"
#
###===  append the thereshold per year

### first split full dataset in years
#cdo splityear ${savepath}/post_cmip6_z500_day.nc ${savepath}/cmip6_z500_day_
#
### now add
#for yy in {1979..2014}
#do
#	cdo merge ${savepath}/cmip6_z500_day_${yy}.nc ${factory}/cmip6_threshold_z500_filtered.nc ${factory}/cmip6_z500_threshold_${yy}.nc
#done

### IDENTIFY REGIONS OF GH ABOVE THE THERESHOLD

#ls ${factory}/cmip6_z500_threshold_* > cmip6_DB_files.txt
#cp cmip6_DB_files.txt cmip6_blocktag_files.txt
#sed -i 's/threshold/blocktag/g' cmip6_blocktag_files.txt
#
#DetectBlobs \
#        --in_data_list cmip6_DB_files.txt \
#        --out_list cmip6_blocktag_files.txt \
#        --thresholdcmd "_DIFF(z,threshold_z),>=,0,0" \
#        --minabslat 25 \
#        --maxabslat 75 \
#        --geofiltercmd "area,>,1e6km2" \
#        --tagvar "block_tag"
#
#### Minimum duration
#cp cmip6_blocktag_files.txt cmip6_blockid_files.txt
#sed -i 's/blocktag/blockid/g' cmip6_blockid_files.txt
#
#StitchBlobs \
#        --in_list cmip6_blocktag_files.txt \
#        --out_list cmip6_blockid_files.txt \
#        --var "block_tag" \
#        --mintime "5d" \
#        --min_overlap_prev 20 \
#        --flatten
#
#### Climatology
Climatology \
        --in_data_list cmip6_blockid_files.txt \
        --out_data ${factory}/cmip6_blocking_climo.nc \
        --var "object_id" \
        --period "seasonal"

