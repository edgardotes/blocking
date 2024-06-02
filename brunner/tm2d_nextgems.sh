#!/bin/bash
#SBATCH --job-name=run_block_rev
#SBATCH -p shared
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -o run_block.out
#SBATCH -e run_block.err
#SBATCH --mem=30Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
#module load python3
module load cdo netcdf-c nco emacs
module load python3/2022.01-gcc-11.2.0

#######################################
###Run ABS blocking detection ICON
######################################
##==extract z500 and compute mean day
#loadpath=/work/bm1235/b382006/regridding/icon/ngc2013_10km_to_1degree/pl/ngc2013_atm_pl_3h*
#savepath=/scratch/b/b382006/tmp/ngc2013/daily

#### create folder to save data
#if [ ! -d ${savepath}  ]; then
#   mkdir -p ${savepath};
#fi
#
#for fn in $loadpath; do
#    fn_new=$savepath/$(basename $fn)
#    output=${fn_new}
#    if [ ! -f output ]; then
#        echo $output
#        cdo -P 48 -daymean -select,name=zfull,level=50000 $fn $output
#   fi
#done

#execute main program *.py input output outname
#python tm2d.icon-ngc2013.py ${savepath} ${savepath} ${savepath}

###===ngc2009
#python tm2d.icon.py

###===ngc3018
#python tm2d.icon_c3.py

############################################
###Run abs blocking detection IFS
################################
#python tm2d.ifs_c3.py

#######################################
###Run ABS blocking detection ERA5
######################################3
#==extract z500 and compute mean day

factory=/scratch/b/b382006/obs/era5/factory/z500
#ln -s /work/bm1235/b382006/era5/z500/era5_regio_prl_195* $factory
#ln -s /work/bm1235/b382006/era5/z500/era5_regio_prl_196* $factory
#ln -s /work/bm1235/b382006/era5/z500/era5_regio_prl_197* $factory
#ln -s /work/bm1235/b382006/era5/z500/era5_regio_prl_198* $factory
loadpath=/scratch/b/b382006/obs/era5/factory/z500/era5_regio_prl_19*
savepath=/scratch/b/b382006/obs/era5/daily/z500
outfile=post_era5_z500_day.nc
#echo $loadpath $outfile

#### create folder to save data
#if [ ! -d ${savepath}  ]; then
#   mkdir -p ${savepath};
#fi
#

#cdo -P 48 -b F32 daymean -mergetime $loadpath $savepath/$outfile

#execute main program *.py input output outname
echo "python tm2d.era5.py ${savepath} ${savepath} ${outfile}"
python tm2d.era5.py ${savepath} ${savepath} ${outfile}
