# PHASTA_TOUS #

The PHASTA_TOUS is a branch of PHASTA solver, which is focused on simulating 
two-phase engineering flows. It is based on the PHASTA_NCSU branch but built to use the open source
svLS linear solver rather than the proprietry Acusim svLS---although the option to use Acusim is
available provided a license is available

For more information, start at our
[wiki page](http://bolotnov.ne.ncsu.edu/index.php?n=Main.HomePage)



### How to compile ###

Regardless of the system you need to make git clone the repository, and then make an empty directory called 'build'.

If you are on insight.ne.ncsu.edu, run

    cmake \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DCMAKE_BUILD_TYPE=Release \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=OFF \
    -DPHASTA_BUILD_PHNSPRE=OFF \
    -DACUSOLVE_LIB=/home/jfang3/svn_develop/develop/LIBLES/lib/x86_64_linux-IB/libles.a \
    ../phasta-ncsu/
    make -j4

Note that the code above may not work for insight as its quite old. 


If you are compiling on intertrack you can create a bash script containing the following:

	#!/bin/bash

	rm -rf build/*  # forcefully remove everything in the build directory
	cd build/       # enter the build directory

	ACUSOLVE_LIB="/home/ctai2/LIBLES:$ACUSOLVE_LIB"
	export ACUSOLVE_LIB


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
	make -j8


The environemnt loaded is paricular to intertrack and a new environemnt will be necessary to run the
code on a different HPC.
As an example the env_intetrack.sh script contains the following:

	#!/bin/bash

	# The required licenses and paths and placed here
	#--------------------------------------------------------------------------------------------------
	export LES_LICENSE_SERVER=multiphase.ne.ncsu.edu

	LD_LIBRARY_PATH="/home/ctai2/LIBLES:$LD_LIBRARY_PATH"
	export LD_LIBRARY_PATH

	ACUSOLVE_LIB="/home/ctai2/LIBLES:$ACUSOLVE_LIB"
	export ACUSOLVE_LIB

	export SIM_LICENSE_FILE=/home/iabolotn/Install/simmetrix/License/NorthCarolinaSU.license
	export MESHSIM=/home/iabolotn/Install/simmetrix/latest
	export SIMMODELER_HOME=/home/iabolotn/Install/simmetrix/SimModeler-latest

	export PATH=/home/iabolotn/Install/simmetrix/SimModeler-latest:$PATH

	#---------------------------------------------------------------------------------------------------

	# purge module
	module purge

	# intel one api
	source /home/ctai2/intel/oneapi/setvars.sh
	# If the above is not included cmake wont be able to find the mpi compile options

	# load gnu9 for mpich
	module load gnu9
	# load mpich for mpicc mpifort mpirun
	module load mpich




