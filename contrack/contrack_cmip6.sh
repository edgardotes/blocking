#!/bin/bash
#SBATCH --job-name=run_block
#SBATCH -p shared
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH -o run_block.out
#SBATCH -e run_block.err
#SBATCH --mem=20Gb
# #SBATCH --reservation=nextGEMS
#SBATCH -A bb1153

###Load modules
#module load python3
#module load python3/2022.01-gcc-11.2.0
#module load cdo netcdf-c nco emacs
#source activate contrack_dev
#source activate python-HH-hackathon
#module use /work/k20200/k202134/hsm-tools/outtake/module
#module add python3/python_3.12-flo
module use /work/k20200/k202134/hsm-tools/outtake/module
module load python3/hamburg-hackathon

#############################################################
###Run ANO Z500 blocking detection cmip6 [MPI-ESM1-2-LR]
#############################################################
###compute anomalies
##python blocks_cmip6_work.py '198?' 'MRI-ESM2-0' 'r1i1p1f1'
##python blocks_cmip6_work.py '199?' 'MRI-ESM2-0' 'r1i1p1f1'
##python blocks_cmip6_work.py '200?' 'MRI-ESM2-0' 'r1i1p1f1'

python blocks_cmip6_work.py '1985' 'MIROC6' 'r1i1p1f1'
python blocks_cmip6_work.py '1995' 'MIROC6' 'r1i1p1f1'
python blocks_cmip6_work.py '2005' 'MIROC6' 'r1i1p1f1'

###compute blockigns
#python contrack.cmip6_z500_anom.py
