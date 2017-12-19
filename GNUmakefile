# Chombo library directory
CHOMBO_HOME   := /home/robertm/Projects/mf-chombo/lib
#CHOMBO_HOME   := /home/robertm/Projects/openStreamer/Chombo-3.2/lib

# Make rules
include $(CHOMBO_HOME)/mk/Make.defs

USE_EB=TRUE
USE_MF=TRUE

# Base file containing int main()
ebase := main

# Requires libraries
LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop

# Target
all: all-test

#
base_dir = .
src_dirs = ./src				\
	./src/amr_mesh				\
	./src/cdr_solver			\
	./src/geometry				\
	./src/global				\
	./src/plasma_solver		 	\
	./src/poisson_solver		 	\
	./src/rte_solver		 	\
	./geometries_prebuilt			\

# Define rules to build everything
include $(CHOMBO_HOME)/mk/Make.example
