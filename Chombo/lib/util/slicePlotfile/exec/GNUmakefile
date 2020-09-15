# -*- Mode: Makefile; -*-

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../

ebase = slicePlot

MINDIM = 2
MAXDIM = 4

# this is the local GNUmakefile which contains this example's multidim
# build info -- libraries and source directory information
MULTIDIM_MAKEFILE = GNUmakefile.multidim

# this is the default target, which just invokes the all-multidim target
all: all-multidim

noLibs: all-multidim-nolibs

# this file contains the basic rules used to build multidim codes (using the 
# GNUmakefile.multidim in this directory), including the shell script
# which orchestrates the make process 
include $(CHOMBO_HOME)/mk/Make.multidim.basic

