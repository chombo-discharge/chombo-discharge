## This file defines variables for use on the login nodes of the NERSC  Cori system (Haswell partition)  
##
## NOTE: everything is always in 64bit mode

makefiles+=local/Make.defs.cori.hsw.gnu

CXX=CC
FC=ftn
MPICXX=CC
USE_64=TRUE

CH_CPP=$(CXX) -E -P

RUN = srun -n 2 ./#

# Compiler flags
cxxoptflags += -O3 -ffast-math -mavx2
foptflags += -O3 -ffast-math -mavx2

# Compile with OpenMP
ifeq ($(OPENMPCC),TRUE)
  cxxoptflags += -fopenmp
  foptflags += -fopenmp
endif

XTRALDFLAGS += -Wl,-zmuldefs

USE_HDF=TRUE

HDFLIBFLAGS=   -L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz
HDFMPILIBFLAGS=-L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz
HDFINCFLAGS=   -I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS) 
HDFMPIINCFLAGS=-I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS)

# NERSC defines its FFTW environment in an odd way...
ifeq ($(USE_FFTW),TRUE)
  FFTWDIR = $(FFTW_INC)/..
endif