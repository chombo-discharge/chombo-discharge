# Chombo library directory
CHOMBO_HOME   := /home/marskar/Projects/mf-chombo/lib
STREAMER_HOME := /home/marskar/chombo-streamer

# Make rules
include $(CHOMBO_HOME)/mk/Make.defs

USE_EB=TRUE
USE_MF=TRUE
DIM=2

# Base file containing int main()
ebase := main

# Requires libraries
LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop

# Target
all: all-test

#
base_dir = .
src_dirs = $(STREAMER_HOME)/src					\
	$(STREAMER_HOME)/src/amr_mesh				\
	$(STREAMER_HOME)/src/cdr_solver				\
	$(STREAMER_HOME)/src/elliptic				\
	$(STREAMER_HOME)/src/geometry				\
	$(STREAMER_HOME)/src/global				\
	$(STREAMER_HOME)/src/plasma_solver		 	\
	$(STREAMER_HOME)/src/poisson_solver		 	\
	$(STREAMER_HOME)/src/rte_solver			 	\
	$(STREAMER_HOME)/src/sigma_solver			\
	$(STREAMER_HOME)/cell_taggers				\
	$(STREAMER_HOME)/geometries_prebuilt			\
	$(STREAMER_HOME)/plasma_models	 			\
	$(STREAMER_HOME)/time_steppers	 			\

# Define rules to build everything
include $(CHOMBO_HOME)/mk/Make.example
