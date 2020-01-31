# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

## Note: requires USE_EB=TRUE

# trace the chain of included makefiles
makefiles += example_EBAMRDiffusion_exec

# the base name(s) of the application(s) in this directory
ebase = diffusionDriver

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
INPUT = diffusion.inputs

include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
