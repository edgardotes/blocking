#!/bin/bash
#SBATCH --job-name=run_block
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
module purge
#module load python3/2022.01-gcc-11.2.0
module load python3
#python3/2022.01-gcc-11.2.0                                       
#python3/2023.01-gcc-11.2.0(default)                              
#python3/unstable

#module load cdo netcdf-c nco emacs

#source activate wb_env

### MAIN DATA PATH
#INPUT=/pool/data/CMIP6/data/CMIP
####

#############################################################
###Run ANO Z500 blocking detection in CMIP6 models
#############################################################
#python tm2d.cmip6.py
python tm2d.cmip6_v2.py
