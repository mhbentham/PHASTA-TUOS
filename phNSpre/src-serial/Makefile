SHOWVERSION :=4.2
ISEXECUTABLE := 1
OUTNAME := NSpre


ifndef ARCHOS
   ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif
SIMARCH := $(ARCHOS)
ifeq ($(ARCHOS),sun4_5) 
    SIMARCH := sun4_5_sc5
endif 
ifeq ($(ARCHOS),x86_64_linux-icc)
    SIMARCH := x64_rhel5_gcc41
    MODELERFLAGS := -lpthread
endif
ifeq ($(ARCHOS),x86_linux-icc) 
    SIMARCH := x86_rhel5_gcc41
    MODELERFLAGS := -lpthread 
endif 
ifeq ($(ARCHOS),x86_linux) 
    SIMARCH := x86_rhel5_gcc41
    MODELERFLAGS := -lpthread 
endif

ifeq ($(ARCHOS),x86_cygwin)
    SIMARCH := x86_win_net
#    MODELERFLAGS := -lpthread  
    ETCFLAGS := $(ETCFLAGS)  ws2_32.lib
endif
dirs := .
#Determine Other Dependancies - you may need to modify this
SHAPES := /usr/local/shapes/latest


ifeq ($(ARCHOS),x86_cygwin)   
    PARASOL_DIR:=C:\program files\Parasolid\Kernel\v16.0\base
else
    SHAPES := /usr/local/shapes/latest
endif

ifndef MESHSIM
  MESHSIM := /users/SCOREC/public/meshSim/latest
#  MESHSIM := /users/SCOREC/public/meshSim/6.1-070807/
endif

ifeq ($(FMDB),1)
  DEFS := $(DEFS) -DFMDB
  FMDBPATH := FMDB
  ifndef MODELER
    MODELER := meshmodel
  endif
else
  DEFS := $(DEFS) -DSIM 
  SIM := 1
  ifndef MODELER
    MODELER := parasolid
  endif
endif

ifeq ($(MODELER),discrete)
  ADDLIB = $(ADDSTDLIB)
  MID := Discrete
  MNAME := Discrete
  DEFS := $(DEFS) -DMODELER_DISCRETE
endif
ifeq ($(MODELER),parasolid)
  MID := parasolid
  MNAME := Parasolid
  CFLAGS := $(CFLAGS) $(MTFLAG)
  CXXFLAGS := $(CXXFLAGS) $(MTFLAG)
  LDFLAGS := $(LDFLAGS) $(MTFLAG)
  INCLUDES := $(INCLUDES) 
  DEFS := $(DEFS) -DMODELER_PARASOLID
endif

ifeq ($(FMDB),1)
  ifeq ($(MODELER),meshmodel)
    MNAME := MeshModel
    DEFS := $(DEFS) -DMESHMODEL
    DEPS := $(DEPS) $(FMDBPATH)/SCORECModel
  endif
endif

OUTNAME := $(OUTNAME)-$(MNAME)
DEFS := $(DEFS) -DVERSION=\"$(SHOWVERSION)\"
#$MESHSIMINTERNAL := /users/SCOREC/public/meshSim/Internal


ifeq ($(ARCHOS),x86_cygwin)
  INCLUDES := $(INCLUDES) -I"$(MESHSIM)/include"
  WIN_IOINCLUDE := $(shell cygpath -w $(DEVROOT)/phasta/phastaIO/include )
  INCLUDES := $(INCLUDES) -I"$(WIN_IOINCLUDE)"
else
 ifeq ($(SIM),1)
  INCLUDES := $(INCLUDES) -I$(MESHSIM)/include 
 endif
  INCLUDES := $(INCLUDES) -Iphasta/phastaIO/
endif

OUTNAME :=$(OUTNAME)_Serial

ifeq ($(ARCHOS),x86_cygwin)
  DEVROOT := $(shell cygpath -w -m $(DEVROOT))
  IOLIB_DIR := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  LIBS := -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib 
else
  DEPS := $(DEPS)  phasta/phastaIO 
  LIBS := $(LIBS) phastaIO 
endif


ifeq ($(MODELER),parasolid)
  MESHSIMMODELVERS = -lSimParasolid220 
  INCLUDES :=$(INCLUDES)
  PARASOLIDPATH := $(MESHSIM)/lib/$(SIMARCH)/psKrnl 
  PSKRNL:= -lpskernel		
  ifeq ($(ARCHOS), sun4_5)
    ETCFLAGS := $(ETCFLAGS) -lF77 -lsunmath  -lsocket -lnsl -lrt
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
     ETCFLAGS := $(ETCFLAGS)   -lnsl   -lpthread 
  endif
  ifeq ($(ARCHOS),x86_linux)
     ETCFLAGS := $(ETCFLAGS)   -lnsl  -lpthread 
  endif
  ifeq ($(ARCHOS),x86_cygwin)
     MESHSIMMODELVERS :=-LIBPATH:"$(MESHSIM)/lib/x86_win_net" SimParasolid160.lib	
     MODELERFLAGS := -LIBPATH:"$(PARASOL_DIR)"  PSKERNEL_ARCHIVE.lib
  endif
endif

ifeq ($(MODELER),discrete)
 ifeq ($(ARCHOS),sgi_6)
   MESHSIMMODELVERS := -lSimDiscrete
 endif
 ifeq ($(ARCHOS),sgi_6-64)
   MESHSIMMODELVERS := -lSimDiscrete
 endif
 ifeq ($(ARCHOS),sun4_5)
   MESHSIMMODELVERS := -lSimDiscrete -lnsl 
   ETCFLAGS := $(ETCFLAGS) -lF77 -lsunmath  -lsocket  -lthread -lrt
 endif
 ifeq ($(ARCHOS),x86_64_linux-icc)
    MESHSIMMODELVERS := -lSimDiscrete
    ETCFLAGS := $(ETCFLAGS)   -lnsl      -lpthread
  endif
 ifeq ($(ARCHOS),x86_linux-icc)
    MESHSIMMODELVERS := -lSimDiscrete
    ETCFLAGS := $(ETCFLAGS)   -lnsl      -lpthread
  endif
 ifeq ($(ARCHOS),x86_linux)
    MESHSIMMODELVERS :=  -lnsl  -lSimDiscrete
    ETCFLAGS := $(ETCFLAGS)      -lpthread
  endif
 ifeq ($(ARCHOS),x86_cygwin)
   MESHSIMMODELVERS := -LIBPATH:"$(MESHSIM)/lib/x86_win_net" SimDiscrete.lib
 endif
endif  
PQUAL :=

ifeq ($(MODELER),meshmodel)
   ifeq ($(ARCHOS),sun4_5)
     ETCFLAGS := $(ETCFLAGS) -lnsl -lF77 -lsunmath -lsocket -lthread -lm
   endif
   ifeq ($(ARCHOS),x86_linux)
      ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread -lm
    endif
   ifeq ($(ARCHOS),x86_linux-icc)
      ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread -lm
    endif
endif

ifeq ($(ARCHOS),x86_cygwin)
  MESHSIMFLAGS :=   $(MESHSIMMODELVERS) \
                    -LIBPATH:"$(MESHSIM)/lib/$(SIMARCH)" SimAdvMeshing.lib \
                     SimModel.lib SimMeshing.lib SimMeshTools.lib \
                     SimModel.lib SimMeshing.lib SimMeshTools.lib
else
   ifeq ($(SIM),1)
   MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH) 
   MESHSIMFLAGS := $(MESHSIMFLAGS)   \
               -lSimAdvMeshing   -lSimMeshing \
                $(MESHSIMMODELVERS)\
                -lSimMeshTools -lSimModel\
		-L$(PARASOLIDPATH) $(PSKRNL)  -lm
   endif
endif

ifeq ($(FMDB),1)
    DEPS := $(DEPS) $(FMDBPATH)/FMDB $(FMDBPATH)/SCORECUtil $(FMDBPATH)/SCORECModel
    LIBS := $(LIBS) FMDB SCORECUtil SCORECModel
endif

ifeq ($(ARCHOS),x86_cygwin)
LDFLAGS := $(LDFLAGS) $(ETCFLAGS) $(MESHSIMFLAGS) $(MODELERFLAGS) $(LIBS)
else
LDFLAGS := $(LDFLAGS) $(ETCFLAGS) $(MESHSIMFLAGS) $(MODELERFLAGS)
LIBS := $(LIBS)  phastaIO
endif


ADDLIB = $(ADDSTDLIB)


# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common

rclean:
	rm .*.d

mclean:
	rm obj/$(ARCHOS)$(QUAL)/mdb2phasta.o
