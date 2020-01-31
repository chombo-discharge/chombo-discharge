# -*- Mode: Makefile;  -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += fourthOrderHyp_advectEqn

# the base name(s) of the application(s) in this directory
ebase = advect

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRTools  BoxTools

# relative paths to source code directories
base_dir = .

FOURTHORDERMAPPED_DIR = $(CHOMBO_HOME)/../example/fourthOrderMappedGrids

src_dirs = ../srcMappedAdvection $(FOURTHORDERMAPPED_DIR)/src $(FOURTHORDERMAPPED_DIR)/util $(FOURTHORDERMAPPED_DIR)/advection $(FOURTHORDERMAPPED_DIR)/advection/mappedCTU $(FOURTHORDERMAPPED_DIR)/example_geometries 

# input file for 'run' target
INPUT = gaussian.inputs.small

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables
#XTRACPPFLAGS += -I$(PACT)/intel/include
#XTRALIBFLAGS +=  -L$(PACT)/intel/lib -lpml -lscore -L../../../ARPACK -larpack_Linux

# application-specific targets
