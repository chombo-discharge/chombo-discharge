# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_MFTools

ebase := testMFStencil testMFpoissonOp testMFVectorPoissonOp testRelaxSolver testBCGMF testBCGVector testMFPoissonOpB

LibNames := MFElliptic EBAMRElliptic AMRElliptic MFTools EBAMRTimeDependent EBAMRTools Workshop EBTools AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ./src

include $(CHOMBO_HOME)/mk/Make.test

