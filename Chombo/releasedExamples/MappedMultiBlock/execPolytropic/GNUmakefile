# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_MappedMultiBlock_execPolytropic

# the base name(s) of the application(s) in this directory
ebase = conservation

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = HOMappedConservation HOMappedGeometry HOAMRTools MOLAMRTimeDependent AMRTimeDependent AMRTools BoxTools BaseTools

# relative paths to source code directories
base_dir = .

POLYTROPIC_DIR = ../../AMRGodunov/srcPolytropic

CARTESIANPOLYTROPIC_DIR = ../../AMRMOL/srcPolytropic

src_dirs = ../src ../srcPolytropic $(CARTESIANPOLYTROPIC_DIR) $(POLYTROPIC_DIR)

# input file for 'run' target
INPUT = gaussian.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# additional libraries

# application-specific variables

# application-specific targets
