## This file defines variables for use on the login nodes of the NERSC Linux
## machine 'hopper'.  
##
## NOTE: everything is always in 64bit mode

makefiles+=local/Make.defs.hopper.gnu

CXX=CC
FC=ftn
MPICXX=CC
#MPICXX=CC -target=linux
USE_64=TRUE

CH_CPP=$(CXX) -march=barcelona -E -P -C

RUN = aprun -n 2 ./#

cxxoptflags += -march=barcelona -ffast-math -O3
foptflags += -O2
# the pgf libs are needed for linking parallel HDF5
flibflags += -lgfortran -L/opt/pgi/default/linux86-64/default/lib 
XTRALDFLAGS += -Wl,-zmuldefs

# The appropriate 'module' must be loaded for this to work.
## dtg: relevant stuff from my .cshrc.ext
## module swap PrgEnv-pgi PrgEnv-gnu
## module load hdf5-parallel
## setenv USE_EB TRUE

USE_HDF=TRUE
#HDF5_DIR = $(CRAY_HDF5_DIR)/hdf5-parallel-gnu

HDFLIBFLAGS=   -L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz -DH5_USE_16_API
HDFMPILIBFLAGS=-L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz -DH5_USE_16_API
HDFINCFLAGS=   -I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS)  -DH5_USE_16_API
HDFMPIINCFLAGS=-I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS)  -DH5_USE_16_API
##
##ifeq ($(USE_64),FALSE)
##  $(error Are you sure you want to run non-64bit?)
##endif
