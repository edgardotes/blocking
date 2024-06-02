#!/bin/bash
#SBATCH --job-name=run_block
#SBATCH -p shared
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -o run_block.out
#SBATCH -e run_block.err
#SBATCH --mem=40Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
module load python3
#module load python3/2022.01-gcc-11.2.0
#module load cdo netcdf-c nco emacs
source activate contrack_dev

#############################################################
###Run ANO Z500 blocking detection ERA
#############################################################
#python contrack.era5_z500.py

###compute anomalies
#python blocks_era5_work.py
###compute blockigns
python contrack.era5_z500_anom.py

#############################################################
###Run ANO VAPV blocking detection ERA
#############################################################
#python contrack.era5_vapv.py
