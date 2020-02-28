#!/bin/csh -fx
# ------------------------------
#SBATCH -A s0818
#SBATCH --export=NONE
#
#PBS -N run
#PBS -o run.log.o%j
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --constraint=hasw
##SBATCH --qos=debug
##PBS -l walltime=1:00:00
##SBATCH --partition=preops
##SBATCH --qos=dastest
##SBATCH --qos=obsdev
##SBATCH --qos=debug
#PBS -l walltime=9:00:00
##PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -j eo
#BSUB -J m2m1c
#BSUB -n 384
#BSUB -W 5:00
#BSUB -o m2m1c.log.o%J
#BSUB -e m2m1c.log.o%J
# ------------------------------
source  /home/bkarpowi/bin/go3.csh
which python3 
./wrfSim_flat_B2B3.py --in $PWD/3km_full_hour_split/wrfout_d01_2019-02-13_00_00_00_section_01 --output-dir $PWD --nthreads 28 --idx 0
