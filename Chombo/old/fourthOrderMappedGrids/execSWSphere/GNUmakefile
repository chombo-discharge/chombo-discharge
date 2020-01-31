# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += fourthOrderHyp_SWEqn

# the base name(s) of the application(s) in this directory
ebase =  swsphere

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRTools  BoxTools

# relative paths to source code directories
base_dir = .

FOURTHORDERMAPPED_DIR = $(CHOMBO_HOME)/../example/fourthOrderMappedGrids

FOURTHORDERHYPERBOLIC_DIR = $(CHOMBO_HOME)/../example/fourthOrderHyperbolic

src_dirs = $(FOURTHORDERMAPPED_DIR)/srcSWSphere $(FOURTHORDERMAPPED_DIR)/srcHyperbolic $(FOURTHORDERMAPPED_DIR)/srcMultiBlock $(FOURTHORDERMAPPED_DIR)/src $(FOURTHORDERMAPPED_DIR)/util $(FOURTHORDERHYPERBOLIC_DIR)/srcCartesian $(CHOMBO_HOME)/../example/InReleasedExamples/AMRGodunov/srcPolytropic

# input file for 'run' target
INPUT = gaussian.inputs.small

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

LIBFLAGS += $(LAPACKLIBS)

# application-specific variables

# application-specific targets
