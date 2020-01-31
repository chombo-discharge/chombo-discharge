# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = diagnostic

## location of 'lib' directory
CHOMBO_HOME = ../../../lib

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames =  AMRTools BoxTools


# input file for 'run' target
INPUT = inputs.compare
# override VERBOSE = # #program doesnt handle -q option

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# include $(CHOMBO_HOME)/mk/Make.test
