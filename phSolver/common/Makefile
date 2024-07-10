SUBSYSNAME := phasta/phSolver/common
MODULENAME := common
NVERS := 1.0
BUILDV := 1.0
NOSHARED = 1
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch )
endif

ifeq ($(ARCHOS),x86_64_linux-icc)
#    INCLUDES := $(INCLUDES) -I/usr/local/apps/mpich2/int101x64/1.3a2/include
#    INCLUDES := $(INCLUDES) -I/usr/local/intel/mpich-1.2.7p1/linux86/10.1/include
    DEFS := $(DEFS) -DLINUX
    LDFLAGS := $(LDFLAGS) -lifcore -lifport -lpthread
    FFLAGS := -W0 -132 
    FFLAGS := $(FFLAGS) -cm  
endif
ifeq ($(ARCHOS),x86_cygwin)
    DEVROOT := $(shell cygpath -w -m $(DEVROOT))
    INCLUDES := $(INCLUDES) -I "C:\Program Files\MPICH\SDK\include"
    DEFS := $(DEFS) -Dintel
endif
ifeq ($(ARCHOS),sun4_5)
    DEFS := $(DEFS) -Dsun4_5 
endif
ifeq ($(ARCHOS),sgi_6-64)
    DEFS := $(DEFS) -Dsgi_6 -Dsgi
endif
ifeq ($(ARCHOS),ia64_linux)
    DEFS := $(DEFS) -Dsun4_5
endif
ifeq ($(ARCHOS),x86_linux-icc)
    DEFS := $(DEFS) -Dsun4_5
    FFLAGS := -w95 -W0 -cm -80
endif
ifeq ($(ARCHOS),alpha_osf)
    FFLAGS := -convert big_endian  -arch host -tune host -unroll 4 -fast -shared -w
    CCFLAGS := -convert big_endian -arch host -tune host -call_shared
    LDFLAGS := -ldxml -lsciport -lc -lfor -lFutil -lm -lots -lUfor
endif
ifeq ($(ARCHOS),ibm_5)
    DEFS := $(DEFS) -Dibm6000
endif
ifeq ($(VERS),)
    DEFS := $(DEFS) -DLAUNCH_GDB
endif

dirs := .
DEPS := $(DEPS) phasta/phastaIO 
mods  := spebc.f pointer.f \
        turbsa.f \
        perprep.f local_mass.f \
        filtprep.f  pvsqbi.f \
        dtn.f bardmc.f \
        timedata.f \
        turbke.f\
        bctint.f readnblk.f spat_var_eps.f\
        bubArr.f

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
