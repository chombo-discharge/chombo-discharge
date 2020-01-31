# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_RichardsonAMRPoisson

# the base name(s) of the application(s) in this directory
ebase = richardsonExtrap

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = OldAMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs =

# input file for 'run' target
#[NOTE: cmd line args do not obey Chombo conventions. <dbs>]
INPUT = -32 inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
