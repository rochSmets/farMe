#!/bin/bash

## ___ sbatch directives
#SBATCH --job-name=weibel
#SBATCH --output=%x-%j.log
#
#SBATCH --ntasks=40
#SBATCH --time=2:00:00
#SBATCH --partition=cpu_shared
#SBATCH --account=phare
##SBATCH --mail-type=ALL
##SBATCH --mail-user=roch.smets@polytechnique.edu

## ___ load modules
module load cmake/3.19.7
module load gcc/10.2.0
module load openmpi/4.1.0
module load hdf5/1.12.1
module load mambaforge/22.11.1-4

#use conda for phare env... created with mamba !
conda activate phare

## ___ set the PYTHONPATH
export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare

## ___ run phare with python
mpirun -n $SLURM_NTASKS python weibel.py

