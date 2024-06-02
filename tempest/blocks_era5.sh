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
loadpath=/work/bm1235/b382006/era5/z500/era5_regio_prl_*
savepath=/scratch/b/b382006/obs/era5/daily/z500
factory=/scratch/b/b382006/obs/era5/factory/z500/tempest

#### create folder to save data
if [ ! -d ${savepath}  ]; then
   mkdir -p ${savepath};
fi

cdo -b F32 daymean -mergetime $loadpath $savepath/era5_z500_day.nc

#python rename_dims_v2.py

###===start running tempest routines
source activate tempest_env

### Compute z500 climatology

###create lists
#ls ${savepath}/post_* > era5_z500_list.txt

#Climatology \
#        --in_data_list era5_z500_list.txt \
#        --out_data ${factory}/era5_z500_LTDM.nc \
#        --var "z" \
#        --period "daily" \
#        --type "mean" \
#        --missingdata
#
#Climatology \
#        --in_data_list era5_z500_list.txt \
#        --out_data ${factory}/era5_z500-2_LTDM.nc \
#        --var "z" \
#        --period "daily" \
#        --type "meansq" \
#        --missingdata

### calculating standard deviation
#
#VariableProcessor \
#        --in_data "${factory}/era5_z500_LTDM.nc;${factory}/era5_z500-2_LTDM.nc" \
#        --out_data "${factory}/era5_z500_mean_stddev.nc" \
#        --var "dailymean_z,_SQRT(_DIFF(dailymeansq_z,_POW(dailymean_z,2)))" \
#        --varout "dailymean_z,stddev_z"
#
#### 4-mode  fourier filter
#FourierFilter \
#        --in_data ${factory}/era5_z500_mean_stddev.nc \
#        --out_data ${factory}/era5_z500_mean_stddev_timesmoothed.nc \
#        --var "dailymean_z,stddev_z" \
#        --dim "time" \
#        --modes 4
#
#FourierFilter \
#        --in_data ${factory}/era5_z500_mean_stddev_timesmoothed.nc \
#        --out_data ${factory}/era5_z500_mean_stddev_smoothed.nc \
#        --var "stddev_z" \
#        --preserve "dailymean_z" \
#        --dim "lon" --modes 2

### computing thereshold
#VariableProcessor \
#        --in_data ${factory}/era5_z500_mean_stddev_smoothed.nc \
#        --out_data ${factory}/era5_threshold_z500_filtered.nc \
#        --var "_SUM(dailymean_z,_MAX(100.0,_PROD(1.5,stddev_z)))" \
#        --varout "threshold_z"

###===  append the thereshold per year

## first split full dataset in years
#cdo splityear post_era5_z500_day.nc era5_z500_day_

##now add
#for yy in {1990..2022}
#do
#cdo merge ${savepath}/era5_z500_day_${yy}.nc ${factory}/era5_threshold_z500_filtered.nc ${factory}/era5_z500_threshold_${yy}.nc
#done

### IDENTIFY REGIONS OF GH ABOVE THE THERESHOLD

#ls ${factory}/era5_z500_threshold_* > era5_DB_files.txt
#cp era5_DB_files.txt era5_blocktag_files.txt
#sed -i 's/threshold/blocktag/g' era5_blocktag_files.txt
#
#DetectBlobs \
#        --in_data_list era5_DB_files.txt \
#        --out_list era5_blocktag_files.txt \
#        --thresholdcmd "_DIFF(z,threshold_z),>=,0,0" \
#        --minabslat 25 \
#        --maxabslat 75 \
#        --geofiltercmd "area,>,1e6km2" \
#        --tagvar "block_tag"
#
### Minimum duration
#cp era5_blocktag_files.txt era5_blockid_files.txt
#sed -i 's/blocktag/blockid/g' era5_blockid_files.txt
#
#StitchBlobs \
#        --in_list era5_blocktag_files.txt \
#        --out_list era5_blockid_files.txt \
#        --var "block_tag" \
#        --mintime "5d" \
#        --min_overlap_prev 20 \
#        --flatten

### Climatology
#Climatology \
#        --in_data_list era5_blockid_files.txt \
#        --out_data ${factory}/era5_blocking_climo.nc \
#        --var "object_id" \
#        --period "seasonal"
#
