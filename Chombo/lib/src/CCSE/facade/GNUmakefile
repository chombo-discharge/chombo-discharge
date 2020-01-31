# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../../..

lbase := facade

LibNames := BoxTools

include $(CHOMBO_HOME)/mk/Make.lib

CCSECPPFLAGS := $(subst -DCH,-DBL,$(CPPFLAGS))
CCSECPPFLAGS := -I.. $(CCSECPPFLAGS)

ofile = o/$(config)/ccseFacade.o

$(ofile) : ccseFacade.cpp
	$(QUIET)$(cxxc) $(CCSECPPFLAGS) $(CXXFLAGS) $(XTRACXXFLAGS) $(CPPFLAGS) $(XTRACPPFLAGS) -DCH_LANG_CC $(cxxcompflag) $< $(cxxobjflag)$@
