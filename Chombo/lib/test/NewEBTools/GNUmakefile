# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

ebase = testIndexedMoments slabTest vofIteratorTest fabCopyTest fabIndexTest ldfabCopyTest fabIOTest testEBAlias EBNormalizeByVolumeFractionTest

LibNames = NewEBAMRTools  NewWorkshop NewEBTools AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test
