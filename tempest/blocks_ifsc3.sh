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

###=== clean tool folder
rm *.txt

##==nextgem model
cycle=3
model=ifs
version=IFS_28-NEMO_25

#==extract z500 and compute mean day
loadpath=/scratch/b/b382006/nextgems/cycle$cycle/$version/z500/z500_${model}_c${cycle}*
savepath=/scratch/b/b382006/nextgems/cycle$cycle/$version/z500/daily
factory=/scratch/b/b382006/nextgems/cycle$cycle/$version/factory/tempest

#### create folder to save data
if [ ! -d ${savepath}  ]; then
   mkdir -p ${savepath};
fi

#cdo -b F32 daymean $loadpath $savepath/${version}_z500_day.nc

###=== prepare for main program
python rename_dims_nextgems.py

###===start running tempest routines
source activate tempest_env

###==== Compute z500 climatology

###create lists
ls ${savepath}/post_* > ${model}_z500_list.txt

Climatology \
        --in_data_list ${model}_z500_list.txt \
        --out_data ${factory}/${model}_z500_LTDM.nc \
        --var "z" \
        --period "daily" \
        --type "mean" \
        --missingdata

Climatology \
        --in_data_list ${model}_z500_list.txt \
        --out_data ${factory}/${model}_z500-2_LTDM.nc \
        --var "z" \
        --period "daily" \
        --type "meansq" \
        --missingdata

###== calculating standard deviation
#
VariableProcessor \
        --in_data "${factory}/${model}_z500_LTDM.nc;${factory}/${model}_z500-2_LTDM.nc" \
        --out_data "${factory}/${model}_z500_mean_stddev.nc" \
        --var "dailymean_z,_SQRT(_DIFF(dailymeansq_z,_POW(dailymean_z,2)))" \
        --varout "dailymean_z,stddev_z"

####=== 4-mode  fourier filter
FourierFilter \
        --in_data ${factory}/${model}_z500_mean_stddev.nc \
        --out_data ${factory}/${model}_z500_mean_stddev_timesmoothed.nc \
        --var "dailymean_z,stddev_z" \
        --dim "time" \
        --modes 4

FourierFilter \
        --in_data ${factory}/${model}_z500_mean_stddev_timesmoothed.nc \
        --out_data ${factory}/${model}_z500_mean_stddev_smoothed.nc \
        --var "stddev_z" \
        --preserve "dailymean_z" \
       --dim "lon" --modes 2

###=== computing thereshold
VariableProcessor \
        --in_data ${factory}/${model}_z500_mean_stddev_smoothed.nc \
        --out_data ${factory}/${model}_threshold_z500_filtered.nc \
        --var "_SUM(dailymean_z,_MAX(100.0,_PROD(1.5,stddev_z)))" \
        --varout "threshold_z"

##########################################
###===  Append the thereshold per year
##########################################

### first split full dataset in years
cdo splityear ${savepath}/post_${version}_z500_day.nc ${savepath}/${version}_z500_day_

#### now add (change years)
for yy in {2020..2025}
do
	cdo merge ${savepath}/${version}_z500_day_${yy}.nc ${factory}/${model}_threshold_z500_filtered.nc ${factory}/${model}_z500_threshold_${yy}.nc
done

####################################################
### IDENTIFY REGIONS OF GH ABOVE THE THERESHOLD
####################################################

ls ${factory}/${model}_z500_threshold_* > ${model}_DB_files.txt
cp ${model}_DB_files.txt ${model}_blocktag_files.txt
sed -i 's/threshold/blocktag/g' ${model}_blocktag_files.txt

DetectBlobs \
        --in_data_list ${model}_DB_files.txt \
        --out_list ${model}_blocktag_files.txt \
        --thresholdcmd "_DIFF(z,threshold_z),>=,0,0" \
        --minabslat 25 \
        --maxabslat 75 \
        --geofiltercmd "area,>,1e6km2" \
        --tagvar "block_tag"

### Minimum duration
cp ${model}_blocktag_files.txt ${model}_blockid_files.txt
sed -i 's/blocktag/blockid/g' ${model}_blockid_files.txt

StitchBlobs \
        --in_list ${model}_blocktag_files.txt \
        --out_list ${model}_blockid_files.txt \
        --var "block_tag" \
        --mintime "5d" \
        --min_overlap_prev 20 \
        --flatten

#### Climatology
Climatology \
        --in_data_list ${model}_blockid_files.txt \
        --out_data ${factory}/${model}_blocking_climo.nc \
        --var "object_id" \
        --period "seasonal"

