# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = loadBalancer

# location of 'lib' directory
CHOMBO_HOME = ../..

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = AMRTools BoxTools

# input file for 'run' target
INPUT = inputs.$(DIM)d.2Proc
override VERBOSE = # #program doesnt handle -q option

include $(CHOMBO_HOME)/mk/Make.test
