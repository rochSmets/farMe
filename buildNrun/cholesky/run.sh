#!/bin/bash

## ___ sbatch directives
#SBATCH --job-name=phare
#SBATCH --output=%x-%j.log
#
#SBATCH --ntasks=200
#SBATCH --time=20:00:00
#SBATCH --partition=cpu_dist
#SBATCH --account=phare
##SBATCH --mail-type=ALL
##SBATCH --mail-user=roch.smets@polytechnique.edu

## ___ load modules
module load cmake/3.19.7
module load gcc/10.2.0
module load openmpi/4.1.0
module load hdf5/1.12.1
module load anaconda3/2020.11

## ___ set the PYTHONPATH
export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare

## ___ run phare with python
mpirun -n  $SLURM_NTASKS  python phare.py

cp alexis.sh ./alexis/
cp alexis.py ./alexis/
