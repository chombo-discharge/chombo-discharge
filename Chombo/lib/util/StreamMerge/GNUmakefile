# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = streamMerge

## location of 'lib' directory
CHOMBO_HOME = ../..

%.template.H : %.template
	xxd -i $< $@
##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = 

override VERBOSE = # #program doesnt handle -q option

include $(CHOMBO_HOME)/mk/Make.test
