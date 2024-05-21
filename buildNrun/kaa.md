
# 15/5/2024

## Preamble

Kaa can be accessed by `ssh`. But the system is not very well updated (I think for reasons related to the [[jupyterhub]]). As a result, gcc version is old, etc.

## Getting modern gcc

By default gcc version is ~8. To get a more recent gcc:

`source /opt/rh/gcc-toolset-11/enable`

Note: there is `gcc-toolset-13` which gives gcc-13 but lacks `gfortran` for some reason.

## Python

There are several version on the system:

- python3.8
- python3.9
- python3.11

### Python virtualenv

I have a [[virtualenv]]:  `~/phareenv2` it uses `python 3.9.17`

### Python packages required

Be sure to have all required packages installed:

```bash
python3 -m pip install -r ~/PHARE/requirements.txt
```


## loading MPI

MPI is not accessible on Kaa unless one does:

```bash
module load mpi
```


## Configuration


Usually when I configure PHARE I do:

```bash
cmake -DCMAKE_CXX_FLAGS="-g3 -O2 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DwithCaliper=OFF -DtestMPI=OFF -DdevMode=OFF
```

But on kaa that does not work and does not find the Parallel version of [[HDF5]]. On kaa there are two hdf5 version of the libraries and headers:

- the sequential version 
	- includes are in `/usr/include` (files: `/usr/include/H5*`)
	- libs are in `/usr/lib64`
- the parallel version:
	- includes are in `/usr/include/openmpi-x86_64
	- libs are in `/usr/lib64/openmpi/lib`

For some reason the `FindHDF5` in PHARE CMake's configuration does not find the parallel version, so we have to provide some defines to CMake. These are :


- `HDF5_INCLUDE_DIRS=/usr/include/openmpi-x86_64`
- `HDF5_C_LIBRARIES=/usr/lib64/openmpi/lib/libhdf5.so `
- `HDF5_LIBRARY_PATH=/usr/lib64/openmpi/lib`

Why the `HDF5_C_LIBRARIES` and not also the `HDF5_CXX_LIBRARIES`? Because in PHARE `HDF5`is required by `HighFive` which is a C++ wrapper over the C HDF5 API and not the C++ one.


So the complete CMake command is (assuming phare source dir is `../pharetest`): 

```bash
cmake -DCMAKE_CXX_FLAGS="-g3 -O2 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -Dphare_configurator=OFF -DwithCaliper=OFF -DtestMPI=OFF -DdevMode=OFF -DHDF5_C_LIBRARIES=/usr/lib64/openmpi/lib/libhdf5.so -DHDF5_INCLUDE_DIRS=/usr/include/openmpi-x86_64 -DHDF5_LIBRARY_PATH=/usr/lib64/openmpi/lib ../pharetest
```


## Building

make -j

