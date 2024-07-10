SUBSYSNAME := phasta/phPost
MODULENAME := Reduce
ISEXECUTABLE := 1
OUTNAME := reduce
dirs := src
DEPS := $(DEPS) phasta/phastaIO
LIBS := phastaIO

ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifeq ($(ARCHOS),x86_cygwin)
  DEVROOT := $(shell cygpath -w -m $(DEVROOT))
  IOLIB_DIR := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  ADDSTDLIB := -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib 
endif

ADDLIB = $(ADDSTDLIB)

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
