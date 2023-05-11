#!/bin/bash

## BEGIN SBATCH directives
#SBATCH --job-name=run047c
#SBATCH --output=run047c.txt
#
#SBATCH --ntasks=80
#SBATCH --time=24:00:00
#SBATCH --partition=cpu_dist
#SBATCH --account=phare
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicolas.aunai@lpp.polytechnique.fr
## END SBATCH directives

## load modules
module load cmake/3.19.7
module load gcc/10.2.0
module load openmpi/4.1.0
module load hdf5/1.12.1
module load mambaforge/22.11.1-4

#use conda for phare env... created with mamba !
conda activate phare

export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare

## execution
mpirun -n $SLURM_NTASKS python harris.py
