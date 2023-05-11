#!/bin/sh

#SBATCH --comment="heckle"
#SBATCH --job-name="blAckDog-b1"
#SBATCH --account=medium
#SBATCH --ntasks=80
#SBATCH --output="hopper.%J.log"
#SBATCH --error="hopper.%J.err"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=roch.smets@lpp.polytechnique.fr

#module load openmpi/gcc/64/1.8.4
#module load phdf5/1.8.16

module load mvapich2/gcc/64/2.1rc1 cmake/3.5.0 hdf5-mvapich/1.10.5

MYRUN=$HOME/shErpA/blAckDog/run/b1/
MYEXE=$HOME/codeS/hecKle/build-heckle/HECKLE

time mpirun -np $SLURM_NTASKS $MYEXE $MYRUN
