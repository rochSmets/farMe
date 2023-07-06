#!/bin/bash

## ___ sbatch directives
#SBATCH --job-name=ionBeam-04a
#SBATCH --output=%x-%j.log
#
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH --partition=cpu_shared
#SBATCH --account=phare
#SBATCH --mail-type=ALL
#SBATCH --mail-user=roch.smets@lpp.polytechnique.fr

## ___ load modules
module load cmake/3.19.7
module load gcc/10.2.0
module load openmpi/4.1.0
module load hdf5/1.12.1
module load mambaforge/22.11.1-4

#use conda for phare env... created with mamba !
conda activate phare

## ___ set the PYTHONPATH
export PYTHONPATH=/mnt/beegfs/workdir/roch.smets/build-release:/mnt/beegfs/home/LPP/roch.smets/codes/far/PHARE/pyphare

## ___ run phare with python
mpirun -n $SLURM_NTASKS python ionBeam-04a.py

