#!/bin/bash
#SBATCH --job-name=run_vapv
#SBATCH -p shared
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -o run_vapv.out
#SBATCH -e run_vapv.err
#SBATCH --mem=20Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
module load python3
module load cdo netcdf-c nco emacs

###=== We are using CDO command

###== save path
##icon ngc2009
#savepath=/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree
###era5
savepath=/work/bm1235/b382006/era5

###== nearest neighbor input
###icon ngc2009
#loadpath=/work/bm1235/b382006/regridding/icon/ngc2009_5km_to_1degree/pl/ngc2009_atm_pl_6h_inst_*.nc
###era5
loadpath=/scratch/b/b382006/era5/pv/era5_regio_prl_19*.nc

###icon ngc2009
#for fn in $loadpath; do
#    fn_new=$(basename $fn)
#    echo $fn_new
#     python blocks_icon_work.py $fn $savepath $fn_new
#done

###era5
for fn in $loadpath; do
    fn_new=$(basename $fn)
    echo $fn_new
     python blocks_era5_work.py $fn $savepath $fn_new
done

###TEST
#infiles=/work/bm1235/b382006/regridding/icon/from_5km_to_1degree/raw_weighted/ngc2009_atm_pl_6h_inst_20210724T000000Z.nc
#fn_new=$(basename $infiles)
#python blocks_icon_work.py $infiles $savepath $fn_new
