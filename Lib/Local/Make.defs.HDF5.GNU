#begin  -- dont change this line

DIM            = 2
DEBUG          = FALSE
OPT            = HIGH
PRECISION      = DOUBLE
PROFILE        = FALSE
CXX            = g++
FC             = gfortran
MPI            = FALSE
MPICXX         = mpic++
CXXSTD         = 14
USE_MT         = TRUE
OPENMPCC       = FALSE
USE_EB         = TRUE
USE_MF         = TRUE
USE_HDF        = TRUE
HDFINCFLAGS    = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
HDFLIBFLAGS    = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5 -lz
HDFMPIINCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
HDFMPILIBFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -lhdf5 -lz
cxxoptflags    = -O3 -march=native
flibflags      = -O3 -march=native
syslibflags    = -llapack -lblas -ldl -lm -lgfortran

#end  -- dont change this line