#!/bin/bash
#SBATCH --job-name=run_block
#SBATCH -p shared
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -o run_block.out
#SBATCH -e run_block.err
#SBATCH --mem=30Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
#module load python3
module load python3/2022.01-gcc-11.2.0
module load cdo netcdf-c nco emacs

#############################################################
###Run ANO Z500 blocking detection IFS
#############################################################
#python contrack.ifs_z500.py
###cycle 3
##===one step
###python contrack.ifs-c3_z500.py

##===two steps
### compute anomalies
python blocks_ifs_work.py

### main program (from anomalies)
#python contrack.ifs-c3_z500_anom.py

#############################################################
###Run ANO Z500 blocking detection ICON
#############################################################
###==Extract geopotential to make easy to load
#loadpath=/work/bm1235/b382006/regridding/icon/ngc2013_10km_to_1degree/pl/ngc2013_atm_pl_3h_inst_20*
#
#savepath=/scratch/b/b382006/nextgems/cycle2/ngc2013/factory/z500
###savepath=/scratch/b/b382006/tmp/ngc2009
#
####== create folder to save data
#if [ ! -d ${savepath}  ]; then
#   mkdir -p ${savepath};
#fi
#
#for fn in $loadpath; do
#    fn_new=$savepath/$(basename $fn)
#    output=${fn_new}
#    if [ ! -f output ]; then
#        echo $output
#        cdo -P 48 -select,name=zfull,level=50000 $fn $output
#
#   fi
#done
#
### compute anomalies
#python blocks_icon_work.py

### main program (from anomalies)
#python contrack.icon_z500_anom.py

###execute main program *.py input output outname
###version computing anomalies at the same time in periods of 10 years
##python contrack.icon_z500.py ${savepath} ${savepath} ${savepath}

#############################################################
###Run ANO VAPV blocking detection ICON
#############################################################
#python contrack.icon_vapv.py
