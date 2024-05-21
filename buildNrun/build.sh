#!/usr/bin/env bash


# get machine name from first word before special character
machine=$(hostname | cut -d'.' -f1 | cut -d'-' -f1)

export PHARE_DIR=$1
export BUILD_DIR=$2

echo "building $machine from $PHARE_DIR in $BUILD_DIR"

#default hdf5 is found by cmake
CMAKE_HDF5=""


# list of machines
machines="kaa cholesky pharemer"

# check if machine in list
if [[ ! " ${machines[@]} " =~ " ${machine} " ]]; then
    echo "Unknown machine"
    exit 1
fi


if [ ${machine}  = "kaa" ]; then
    source ~/phare_env2/bin/activate
    source /opt/rh/gcc-toolset-11/enable # toolset-13 lacks gfortran
    CMAKE_HDF5="-DHDF5_C_LIBRARIES=/usr/lib64/openmpi/lib/libhdf5.so -DHDF5_C_INCLUDE_DIRS=/usr/include/openmpi-x86_64  -DHDF5_LIBRARY_PATH=/usr/lib64/openmpi/lib/"
fi

cd $BUILD_DIR
echo $PHARE_DIR
CMAKE_CXX_FLAGS="-DNDEBUG -g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1"
cmake  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" -Dphare_configurator=OFF -DtestMPI=OFF -DdevMode=OFF  ${CMAKE_HDF5} -DwithCcache=0 ${PHARE_DIR} 
make -j


#cd $HOME
#[ ! -d "PHARE" ] && git clone https://github.com/PHAREHUB/PHARE --recursive
#
#cd PHARE
#[ ! -d ".venv" ] && python3 -m venv .venv
#. .venv/bin/activate
#[ ! -f ".pip_installed" ] && python3 -m pip install -r requirements.txt && echo 1 > .pip_installed
#
#export PYTHONPATH="${WORK}/build:${PWD}:${PWD}/pyphare"
#
#cd $WORK
#
## write script it missing
#[ ! -f "build.sh" ] && cat > build.sh << EOL
#mkdir -p build
#cd build
#CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native" # -DPHARE_DIAG_DOUBLES=1 
#cmake ~/PHARE/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="\${CMAKE_CXX_FLAGS}" -Dphare_configurator=ON
#make
#EOL
#chmod +x build.sh
#[ ! -d "build" ] && ./build.sh
#


