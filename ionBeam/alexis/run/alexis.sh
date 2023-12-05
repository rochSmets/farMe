#!/bin/bash

## ___ sbatch directives
#SBATCH --job-name=alexis
#SBATCH --output=%x-%j.log
#
#SBATCH --ntasks=40
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
#module load anaconda3/2020.11
module load python/3.9.2

source $HOME/venv/phare/bin/activate

## ___ set the PYTHONPATH
export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare

## ___ run phare with python
mpirun -n  $SLURM_NTASKS  python alexis.py

