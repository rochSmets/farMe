#!/bin/bash
# ccc_msub -m scratch,store irene.sh
# ccc_mpp -u $USER


RUN=blAckDog/run/07a/

MYRUN=$HOME/shErpA/$RUN
MYEXE=$HOME/codeS/hecKle/heckle/bin/heckle.exe
STORE=$CCCSTOREDIR/$RUN

# ------------------------------------------------------------------------------
#MSUB -q knl            # Partition
#MSUB -A gen10417       # Project ID
#MSUB -r blackDog       # Request name
#MSUB -T 86400          # Elapsed time limit in seconds : 24 hours
#MSUB -n 6400           # Total number of tasks to use
#MSUB -o out.%I.log     # Standard output. %I is the job id
#MSUB -e err.%I.log     # Error output. %I is the job id
#MSUB -@ roch.smets@lpp.polytechnique.fr:begin,end
# ------------------------------------------------------------------------------

set -x

cd $CCCSCRATCHDIR

module sw feature/openmpi/net/bxi
module load flavor/hdf5/parallel
module load hdf5

ccc_mprun $MYEXE $MYRUN

mv $BRIDGE_MSUB_PWD/*.log $MYRUN

mv *.h5 $STORE
