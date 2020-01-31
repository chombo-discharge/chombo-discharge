# -*- Mode: Makefile -*- 

## This makefile builds one or more subdirectories.

## NOTE: requires USE_EB=TRUE

## Define the variables needed by Make.example.subdirs

# trace the chain of included makefiles
makefiles += example_RichardsonAMRPoisson

# relative path to the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# relative paths to the subdirectories to build
exec_dirs = exec

# application-specific variables go here ...

# the shared code for building in example subdirectories
include $(CHOMBO_HOME)/mk/Make.example.subdirs

# application-specific targets go here ...
