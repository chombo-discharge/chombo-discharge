# -*- Mode: Makefile; Modified: "Fri 19 Aug 2005 16:02:15 by dbs"; -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
# makefiles += fourthOrderHyp_advectEqn

# the base name(s) of the application(s) in this directory
ebase = conservation

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# added by petermc, 15 Dec 2008
LAPACKLIBS = -llapack # -lblas

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .

FOURTHORDERMAPPED_DIR = $(CHOMBO_HOME)/../example/fourthOrderMappedGrids

# POLYTROPIC_DIR = $(CHOMBO_HOME)/../example/MappedGridPPMAMRGodunov/srcPolytropic
POLYTROPIC_DIR = $(CHOMBO_HOME)/../example/AMRGodunov/srcPolytropic

src_dirs = ../srcMapped ../srcCartesian $(FOURTHORDERMAPPED_DIR)/src $(FOURTHORDERMAPPED_DIR)/example_geometries $(FOURTHORDERMAPPED_DIR)/util $(POLYTROPIC_DIR)

# input file for 'run' target
INPUT = gaussian.inputs.small

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# added by petermc, 15 Dec 2008
LIBFLAGS += $(LAPACKLIBS)

# application-specific variables

# application-specific targets
