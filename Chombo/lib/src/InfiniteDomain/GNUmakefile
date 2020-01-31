# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := infinitedomain

LibNames := BoxTools

#ifneq ($(USE_FFTW)$(DIM),TRUE3)
#  $(warning warning: InfiniteDomain requires USE_FFTW=TRUE and DIM=3 to build)
#else
#  lbase := dummyinfinitedomain
#endif

include $(CHOMBO_HOME)/mk/Make.lib





