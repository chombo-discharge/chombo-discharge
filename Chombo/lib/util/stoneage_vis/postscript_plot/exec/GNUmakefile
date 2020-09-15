# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_AMRGodunov_execPolytropic

# the base name(s) of the application(s) in this directory
ebase = plotFromHDF

# the location of Chombo lib dir
CHOMBO_HOME = $(HOME)/_svnchombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames =  EBTools  AMRTimeDependent AMRTools  BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
INPUT = ramp.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
