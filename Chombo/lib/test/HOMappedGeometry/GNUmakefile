# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_homappedgeometry

ebase := testOldMBFR testCartesian testRThetaZ cubedSphere2DConst \
	cubedSphere2DTest testMultiBlockFluxRegister \
	testCubedSphereBlockRegister testMBAggStencil

# These take way too long.
# ebase += testMBLevelExchange testMBLevelCopier

LibNames := HOMappedGeometry MOLAMRTimeDependent AMRTimeDependent HOAMRTools AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ./src

include $(CHOMBO_HOME)/mk/Make.test

LIBFLAGS += $(LAPACKLIBS)
