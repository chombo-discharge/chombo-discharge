# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += multiblock_test

# the base name(s) of the application(s) in this directory
ebase = testMBLevelExchange


# the location of Chombo lib dir
CHOMBO_HOME = $(shell pwd)/../../../lib


# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTools BoxTools

# relative paths to source code directories
base_dir = .


src_dirs = ../src ../example_geometries ../srcMultiBlock ../util ../srcTest

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
