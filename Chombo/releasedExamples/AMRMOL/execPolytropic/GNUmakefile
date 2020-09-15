# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_AMRMOL_execPolytropic

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = conservation

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = MOLAMRTimeDependent AMRTimeDependent HOAMRTools AMRTools BoxTools


# relative paths to source code directories
base_dir = .

POLYTROPIC_DIR = ../../AMRGodunov/srcPolytropic

src_dirs = ../srcPolytropic $(POLYTROPIC_DIR)

# input file for 'run' target
INPUT = gaussian.inputs.small

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# additional libraries

# application-specific variables

# application-specific targets
