#!/bin/bash

#rm -rf build/*	# forcefully remove everything in the build directory
cd build/	# enter the build directory

#module load intel/2022.1.0
ACUSOLVE_LIB="/home/ctai2/LIBLES:$ACUSOLVE_LIB"
export ACUSOLVE_LIB

#source /home/ctai2/intel/oneapi/setvars.sh

#Need to source an environment with the correct compilers. This step will
#depend on the compilers you have available on your machine
source /home/mhbentham/development/envs/env_intertrack.sh


cmake \
-DCMAKE_C_COMPILER=icc \
-DCMAKE_CXX_COMPILER=icpc \
-DCMAKE_Fortran_COMPILER=ifort \
-DCMAKE_BUILD_TYPE=Release \
-DPHASTA_INCOMPRESSIBLE=ON \
-DPHASTA_COMPRESSIBLE=OFF \
-DPHASTA_BUILD_PHNSPRE=ON \
-DACUSOLVE_LIB=/home/ctai2/LIBLES/libles.a \
-DPHASTA_SOURCE_DIRECTORY=. \
..
#-B . -S ..
make -j8

