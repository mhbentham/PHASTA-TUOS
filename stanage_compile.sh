#!/bin/bash

rm -rf build/*	# forcefully remove everything in the build directory
cd build/	# enter the build directory

module load impi/2021.7.1-intel-compilers-2022.2.1

# set the location of the acusolve_lib. Note that this library shouldnt actually
# be used as the svLS library is used instead. svLS is open source, whereas
# acusolve is proprietry.

# YOU WILL NEED TO CHANGE THE PATH BELOW
ACUSOLVE_LIB="/users/mep23mhb/PHASTA/misc/stanage_setup/LIBLES:$ACUSOLVE_LIB"
export ACUSOLVE_LIB


#Need to source an environment with the correct compilers. This step will
#depend on the compilers you have available on your machine
#source /home/mhbentham/development/envs/env_intertrack.sh


cmake \
-DCMAKE_C_COMPILER=mpiicc \
-DCMAKE_CXX_COMPILER=mpiicpc \
-DCMAKE_Fortran_COMPILER=mpiifort \
-DCMAKE_BUILD_TYPE=Debug \
-DPHASTA_INCOMPRESSIBLE=ON \
-DPHASTA_COMPRESSIBLE=OFF \
-DPHASTA_BUILD_PHNSPRE=ON \
-DACUSOLVE_LIB=/users/mep23mhb/PHASTA/misc/stanage_setup/LIBLES/libles.a \
-DPHASTA_SOURCE_DIRECTORY=. \
..
#-B . -S ..
make -j8

