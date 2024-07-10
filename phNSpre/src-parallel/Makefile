SUBSYSNAME := phasta
MODULENAME := phNSpre2
BUILDV := 1.0
ISEXECUTABLE := 1
OUTNAME := NSpre_Parallel
F77 := f90

ifndef ARCHOS
      ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

dirs := .

ifeq ($(ARCHOS), x86_linux-icc)
  ADDSTDLIB := -lifcoremt -limf -lifport -lpthread
endif
ifeq ($(ARCHOS), x86_64_linux-icc)
  ADDSTDLIB := -lifcoremt -limf -lifport -lpthread
endif
ifeq ($(ARCHOS), x86_linux)
  ADDSTDLIB := -lg2c
endif
ifeq ($(ARCHOS), sgi_6-64)
  ADDSTDLIB := -lftn -lfortran -lm 
  DEFS := -Dsgi 
endif
ifeq ($(ARCHOS), sgi_6)
  ADDSTDLIB := -lftn -lfortran -lm 
  DEFS := -Dsgi 
endif
ifeq ($(ARCHOS),sun4_5)
  ADDSTDLIB :=-lfsu -lfai2 -lfsumai -lfminlai -lfminvai -lfmaxvai -lfmaxlai -lF77 -lsunmath
  CXXFLAGS := -w
  DEFS := -Dsun4_5 
  MTFLAG := -mt
endif
ifeq ($(ARCHOS),sun4_5-64)
  ADDSTDLIB :=-lfsu -lfai2 -lfsumai -lfminlai -lfminvai -lfmaxvai -lfmaxlai -lF77 -lsunmath
  CXXFLAGS := -w
  DEFS := -Dsun4_5 
  MTFLAG := -mt
endif
     

ADDLIB = $(ADDSTDLIB)


ifeq ($(ARCHOS),x86_cygwin)
  INCLUDES       := $(INCLUDES) -I"$(MESHSIM)/include"
  WIN_IOINCLUDE  := $(shell cygpath -w $(DEVROOT)/phasta/phastaIO/include )
  INCLUDES       := $(INCLUDES) -I"$(WIN_IOINCLUDE)"
  METIS_INC      := $(shell cygpath -w $(DEVROOT)/phasta/phMetis/include )
  INCLUDES       := $(INCLUDES) -I"$(METIS_INC)"
else
  INCLUDES := $(INCLUDES) -I$(MESHSIM)/include 
  INCLUDES := $(INCLUDES) -Iphasta/phastaIO/
  INCLUDES := $(INCLUDES) -Iphasta/phMetis
endif




ifeq ($(ARCHOS),x86_cygwin)
  METISLIB_DIR := $(DEVROOT)/phasta/phMetis/lib/$(ARCHOS)
  WINMETISLIB_DIR := $(shell cygpath -w $(METISLIB_DIR))
  DEPS :=  -LIBPATH:"$(WINMETISLIB_DIR)" libmetis.lib
  IOLIB_DIR := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  DEPS := $(DEPS) -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib
  ADDSTDLIB:=  -LIBPATH:"c:\program files\intel\compiler\fortran\9.0\ia32\lib"
else
DEPS := $(DEPS) \
        phasta/phMetis \
        phasta/phastaIO
endif


ifeq ($(ARCHOS),x86_cygwin)
    ADDLIB := $(DEPS) $(ADDSTDLIB)
endif


LIBS := $(LIBS) metis phastaIO

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common

rclean:
	rm .*.d
