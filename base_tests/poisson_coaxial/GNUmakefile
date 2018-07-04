# Chombo and chombo-streamer directories 
CHOMBO_HOME   := /home/robertm/Projects/mf-chombo/lib
STREAMER_HOME := /home/robertm/Projects/plasmac

include $(CHOMBO_HOME)/mk/Make.defs

USE_EB=TRUE
USE_MF=TRUE
DIM=2


# Base file containing int main()
ebase := main

LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop


# Target
all: all-test

base_dir = .
src_dirs = $(STREAMER_HOME)/src \
	$(STREAMER_HOME)/src/amr_mesh \
	$(STREAMER_HOME)/src/cdr_solver \
	$(STREAMER_HOME)/src/elliptic \
	$(STREAMER_HOME)/src/geometry \
	$(STREAMER_HOME)/src/global \
	$(STREAMER_HOME)/src/plasma_solver \
	$(STREAMER_HOME)/src/poisson_solver \
	$(STREAMER_HOME)/src/rte_solver \
	$(STREAMER_HOME)/src/sigma_solver \
	$(STREAMER_HOME)/geometries_prebuilt/coaxial_cable\
	$(STREAMER_HOME)/time_steppers/rk2\
	$(STREAMER_HOME)/plasma_models/advection_kinetics


include $(CHOMBO_HOME)/mk/Make.example

