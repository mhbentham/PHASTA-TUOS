project(phSolver)
set(CMAKE_Fortran_MODULE_DIRECTORY ${phSolver_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(CASES "/home/ssatpat/SyncIO_test"
  CACHE string
  "/home/ssatpat/SyncIO_test")
option(PHASTA_TESTING "enable tests" ON)

add_subdirectory(common)
option(PHASTA_INCOMPRESSIBLE "Build the incompressible solver" ON)
option(PHASTA_COMPRESSIBLE "Build the compressible solver" OFF)
if(PHASTA_INCOMPRESSIBLE)
add_subdirectory(incompressible)
endif(PHASTA_INCOMPRESSIBLE)
if(PHASTA_COMPRESSIBLE)
add_subdirectory(compressible)
endif(PHASTA_COMPRESSIBLE)
