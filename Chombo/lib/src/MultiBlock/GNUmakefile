# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

# added by petermc, 7 Apr 2008
# LAPACK_HOME = /home/bergman/users/petermc/LAPACK/current
# ifeq ($(DEBUG),TRUE)
#    LAPACKLIBS = $(LAPACK_HOME)/lapack_LINUX_DEBUG.a $(LAPACK_HOME)/blas_LINUX_DEBUG.a $(LAPACK_HOME)/tmglib_LINUX_DEBUG.a
# else
#    LAPACKLIBS = $(LAPACK_HOME)/lapack_LINUX.a $(LAPACK_HOME)/blas_LINUX.a $(LAPACK_HOME)/tmglib_LINUX.a
# endif
# 
# XTRALIBFLAGS += $(LAPACKLIBS)

lbase := multiblock

LibNames := AMRElliptic AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.lib
