# PHASTA_NCSU #

The PHASTA_NCSU is a branch of PHASTA solver, which is focused on simulating 
two-phase engineering flows. 

For more information, start at our
[wiki page](http://bolotnov.ne.ncsu.edu/index.php?n=Main.HomePage)


Below are the old instructions to compile the code for reference.
### How to compile ###

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

If you are on Mira, Cetus or cooley, run

    cmake \
    -DCMAKE_TOOLCHAIN_FILE=/home/fang/CMakeToolchainFiles/BlueGene/Q/BGQ-XLMPI-toolchain.cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_DEBUG="-O2 -qstrict" \
    -DCMAKE_C_FLAGS_DEBUG="-O2 -qstrict" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-O2 -qstrict" \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=OFF \
    -DACUSOLVE_LIB=/home/fang/LIBLES/Mira/libles-xl-dbg-v1r2m2-150330-feb2015cmp.a  \
    ../phasta-ncsu
    make -j8
