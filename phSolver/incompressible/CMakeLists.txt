# 	NOTE: it is likely this file will have to be edited heavily to allow the
# 	svLS version of the code to work as the structure between the NCSU and 
# 	UOCB version are very different. This file will therefore be the cause of
# 	plenty of errors - Magnus

file(GLOB SRC_FILES *.cc *.f *.f90 *.c)
list(REMOVE_ITEM SRC_FILES "main.c")
add_library(incompressible ${SRC_FILES})
find_library(ACUSOLVE_LIB libles)

# BEGIN MAGNUS------------------------------------------------------------------
# Block prepares for using svLS by setting a definition if enabled, and
# sending a message if it has not been enabled
if(PHASTA_USE_SVLS)
	add_definitions(-DHAVE_SVLS)
endif()

if(NOT PHASTA_USE_SVLS)
	message(FATAL_ERROR
	"You have not enabled svLS. You must enable it via \
	PHASTA_USE_SVLS")
endif()

# add the LES_LIB FLAG for the acusim solver
if (PHASTA_USE_LESLIB)
	add_definitions(-DHAVE_LESLIB)
endif()
# END MAGNUS------------------------------------------------------------------

if(PHASTA_USE_MPI)
	include_directories(${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
endif(PHASTA_USE_MPI)
include_directories(${PHASTA_BINARY_DIR}/include ${PHASTA_BINARY_DIR})
add_executable(phastaIC.exe main.c)
set_target_properties(phastaIC.exe PROPERTIES HAS_CXX TRUE)
set_target_properties(incompressible PROPERTIES HAS_CXX TRUE)
set_target_properties(phastaIC.exe PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(incompressible PROPERTIES LINKER_LANGUAGE Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
	#force preprocessing of solfar.f
	set_source_files_properties(solfar.f PROPERTIES COMPILE_FLAGS "-qsuffix=cpp=f")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	#force preprocessing of solfar.f
	set_source_files_properties(solfar.f PROPERTIES COMPILE_FLAGS "-Mpreprocess")
endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#incompressible reuses a fortran module from common w/o telling cmake
#this is why we specify a module directory in the next cmakelists up
#and we need to serialize incompressible after common
add_dependencies(incompressible common)
if(PHASTA_USE_MPI)
	target_link_libraries(phastaIC.exe ${MPI_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) 
endif(PHASTA_USE_MPI)

# BEGIN MAGNUS-----------------------------------------------------------------------
#  	this block is a work in process, some of the variables have not been declared yet
# 	and so will likely lead to errors
# The UOCB links the libraries using the PHASTA_COMMON_LIBS library as such:

#if(PHASTA_USE_SVLS)
#	message(STATUS "using svLS")
#	set(PHASTAIC_LIBS
#	${PHASTA_COMMON_LIBS}
#	incompressible
#	svLS
#	incompressible
#	${PHASTA_COMMON_LIBS})
#endif()
# END MAGNUS--------------------------------------------------------------------------

target_link_libraries(phastaIC.exe ${ACUSOLVE_LIB})

# BEGIN MAGNUS---------------------------------------------------------------------
#However, we havent defined this variable, instead we will try to link the libraries using the following:
target_link_libraries(phastaIC.exe svLS)
#target_link_libraries(phastaIC.exe svLS incompressible)
# END MAGNUS-----------------------------------------------------------------------

target_link_libraries(phastaIC.exe common incompressible phMetis common phMetis phastaIO)
if(PHASTA_BUILD_SHAPEFUNCTION)
	target_link_libraries(phastaIC.exe shapeFunction)
endif(PHASTA_BUILD_SHAPEFUNCTION)
if(PHASTA_BUILD_PHSHAPE)
	target_link_libraries(phastaIC.exe phshape)
endif(PHASTA_BUILD_PHSHAPE)
target_link_libraries(phastaIC.exe ${ACUSOLVE_LIB} incompressible)

#-----------------------------------------------------------------
# BEGIN MAGNUS
# link svLS library
target_link_libraries(phastaIC.exe svLS incompressible)
# END MAGNUS
#-----------------------------------------------------------------

if(PHASTA_TESTING)
  add_subdirectory(test)
endif(PHASTA_TESTING)

