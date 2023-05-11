#!/bin/bash

## ___ sbatch directives
#SBATCH --job-name=compil_far
#SBATCH --output=%x-%j.log
#
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=00:30:00
#SBATCH --partition=cpu_test
#SBATCH --account=phare
##SBATCH --mail-type=ALL
##SBATCH --mail-user=roch.smets@polytechnique.edu

## ___ load modules
# module load cmake/3.19.7
# module load gcc/10.2.0
# module load openmpi/4.1.0
# module load hdf5/1.12.1
# module load anaconda3/2020.11

## ___ set the PYTHONPATH
# export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare

## ___ compile phare
# conda activate phare
# cmake-phare
cd $HOME/codes/far/build-release
make -j40
