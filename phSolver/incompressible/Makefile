SUBSYSNAME := phasta/phSolver/incompressible
MODULENAME := incompressible
BUILDV := 1.0
ISEXECUTABLE := 1
OUTNAME := phastaIC.exe
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

  PRECOMPILE := copymods
  ARPACKPATH := $(DEVROOT)/ARPACK

  LIBARPACK := $(ARPACKPATH)/libarpack_linux.a
  LIBPARPACK := $(ARPACKPATH)/parpack_MPI-linux.a
  
  ifeq ($(LESLIBDIR), )
       LESLIBDIR := /Install/develop/LIBLES
  endif
  ifeq ($(ARCHOS),x86_64_linux-icc)
#	  CXX := mpicxx
	  INCLUDES := $(INCLUDES) -I/usr/local/include
	  DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
	  CFLAGS := -w
	  FFLAGS := -W0 -132 -fpp
#          FFLAGS := -ffixed-line-length-132
	  ifeq ($(AMG),1)
	      LDFLAGS := $(LDFLAGS) $(LIBPARPACK) $(LIBARPACK)
	  endif
#          LDFLAGS := $(LDFLAGS) -L/usr/local/apps/mpich2/int101x64/1.3a2/lib -lmpl -L/usr/local/intel/fce/10.1.022/lib -lpthread -gfortran
#          LDFLAGS := $(LDFLAGS) -L/usr/local/lib -lmpl -L/usr/local/intel/lib -lmpl -lifcore -lifport -lpthread -lmpl
          LDFLAGS := $(LDFLAGS) -lifcore  -lifport -lpthread -lmpl
#          LDFLAGS := $(LDFLAGS) -L/usr/local/apps/mpich2/int101x64/1.3a2/lib -lmpl -L/usr/local/intel/fce/10.1.022/lib -lifcore -lifport -lpthread -L/usr/local/apps/mpich2/int101x64/1.3a2/lib -lmpl -lm
#          LDFLAGS := $(LDFLAGS) -L/usr/local/intel/fce/10.1.022/lib -lifcore -lifport -lpthread -L/usr/local/apps/mpich2/int101x64/1.3a2/lib -lmpl
#	  LDFLAGS := $(LDFLAGS) -lifcore -lifport -lpthread -lmpl
#	  LDFLAGS := $(LDFLAGS) -L/usr/local/intel/mkl72/lib/32 -lmkl_ia32 -lmkl_p4 -lmkl_lapack
#	  LDFLAGS := -lintrins -lCEPCF90 -lPEPCF90 -lF90
#         SIMARCH := amd64_rhel3_gcc32
          SIMARCH := amd64_rhel3_gcc32
  endif

ifeq ($(AMG),1)
    DEFS := $(DEFS) -DAMG
endif


dirs := .
mods  := stats.f
DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver phasta/phastaIO phasta/phMetis\
        phasta/shapeFunction LIBLES/1.5 
LIBS := common phastaIO metis shapeFunction les

ifeq ($(AMG),1)
    DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver
	LIBS := $(LIBS) AMG
endif

INCLUDES := $(INCLUDES) -I../common
ADDLIB = $(ADDSTDLIB)
# need to include common for intel fortran .mod files
include $(DEVROOT)/Util/buildUtil/make.common

.PHONY: copymods
.PHONY: $(OUTFILE)

copymods:
	cp ../common/*.mod .
ifeq ($(AMG),1)	
	cp ../AMG/*.mod .
endif
