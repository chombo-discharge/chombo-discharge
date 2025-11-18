#begin  -- dont change this line

DIM            = 2
DEBUG          = FALSE
OPT            = HIGH
PRECISION      = DOUBLE
PROFILE        = FALSE
CXX            = g++-6
FC             = gfortran
MPI            = TRUE
MPICXX         = mpicxx
CXXSTD         = 14
USE_PETSC      = FALSE
USE_EB         = TRUE
USE_HDF        = TRUE
HDFINCFLAGS    = -I/usr/local/hdf5-1.8.14/include
HDFLIBFLAGS    = -L/usr/local/hdf5-1.8.14/lib -lhdf5 -lz
HDFMPIINCFLAGS = -I/usr/local/hdf5-1.8.14/include
HDFMPILIBFLAGS = -L/usr/local/hdf5-1.8.14/lib -lhdf5 -lz
USE_MF         = TRUE
cxxoptflags    = -O3
flibflags      = -O3
syslibflags    = -llapack -lblas -ldl -lm

#end  -- dont change this line
