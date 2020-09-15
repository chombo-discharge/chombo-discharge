# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_MappedMultiBlock_execAdvection

# the base name(s) of the application(s) in this directory
ebase = advect

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = HOMappedConservation HOMappedGeometry MOLAMRTimeDependent HOAMRTools AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .

src_dirs = ../src ../srcAdvection

# input file for 'run' target
INPUT = cosinebelladvection.cubedsphere2d.064.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# additional libraries

# application-specific variables

# application-specific targets
