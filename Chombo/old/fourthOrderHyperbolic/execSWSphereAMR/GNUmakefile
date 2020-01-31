# -*- Mode: Makefile -*- 

USE_PYTHON=TRUE

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_fourthOrderHyperbolic_execSWSphereAMR

# the base name(s) of the application(s) in this directory
ebase = swsphere

# the location of Chombo lib dir
CHOMBO_HOME = $(shell pwd)/../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = PyParse AMRTimeDependent AMRTools BoxTools BaseTools

# relative paths to source code directories
base_dir = .

FOURTHORDERMAPPED_DIR = $(CHOMBO_HOME)/../example/fourthOrderMappedGrids
POLYTROPIC_DIR = $(CHOMBO_HOME)/../example/AMRGodunov/srcPolytropic

# Note: we have to include the Polytropic directory here because the 
# mapped AMR directory needs the gasdynamic solver.
src_dirs = ../srcMappedAMR $(POLYTROPIC_DIR) $(FOURTHORDERMAPPED_DIR)/src $(FOURTHORDERMAPPED_DIR)/util $(FOURTHORDERMAPPED_DIR)/srcMultiBlock $(FOURTHORDERMAPPED_DIR)/srcMultiBlock/AMRTools

# input file for 'run' target
INPUT = gaussian.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# additional libraries
# LAPACKLIBS = -llapack
LIBFLAGS += $(LAPACKLIBS)

# application-specific variables

# application-specific targets
