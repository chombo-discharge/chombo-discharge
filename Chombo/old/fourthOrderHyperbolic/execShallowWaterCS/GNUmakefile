# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += SWSphere

# the base name(s) of the application(s) in this directory
ebase =  swsphere


# the location of Chombo lib dir
CHOMBO_HOME = $(shell pwd)/../../../lib

LAPACKLIBS = -llapack -lblas

# names of Chombo libraries needed by this program, in order of search.
LibNames = PyParse MultiBlock AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
FOURTHORDERMAPPED_DIR = $(CHOMBO_HOME)/../example/fourthOrderMappedGrids
FOURTHORDERHYPERBOLIC_DIR = $(CHOMBO_HOME)/../example/fourthOrderHyperbolic

src_dirs = ../srcSWSphere ../src \
					$(FOURTHORDERMAPPED_DIR)/srcHyperbolic \
					$(FOURTHORDERMAPPED_DIR)/srcAdvection \
					$(FOURTHORDERHYPERBOLIC_DIR)/srcCartesian \
					$(FOURTHORDERMAPPED_DIR)/srcSWSphere \
					$(FOURTHORDERMAPPED_DIR)/srcMultiBlock \
					$(FOURTHORDERMAPPED_DIR)/srcMultiBlock/AMRTools \
					$(FOURTHORDERMAPPED_DIR)/util \
					$(CHOMBO_HOME)/../example/AMRGodunov/srcPolytropic # Yuck! Needed by 4th-order Cartesian stuff.

# input file for 'run' target
INPUT = inputs


# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

LIBFLAGS += $(LAPACKLIBS)

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
