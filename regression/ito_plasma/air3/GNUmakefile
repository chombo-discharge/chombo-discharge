export CHOMBO_HOME=$(PLASMAC_HOME)/Chombo/lib
include $(CHOMBO_HOME)/mk/Make.defs

USE_EB=TRUE
USE_MF=TRUE
DIM=2


# Base file containing int main()
ebase := main

LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools


# Target
all: all-test

base_dir = .
src_dirs = $(PLASMAC_HOME)/src \
	$(PLASMAC_HOME)/src/amr_mesh \
	$(PLASMAC_HOME)/src/cdr_solver \
	$(PLASMAC_HOME)/src/cell_tagger \
	$(PLASMAC_HOME)/src/elliptic \
	$(PLASMAC_HOME)/src/geometry \
	$(PLASMAC_HOME)/src/global \
	$(PLASMAC_HOME)/src/driver \
	$(PLASMAC_HOME)/src/particle \
	$(PLASMAC_HOME)/src/ito_solver \
	$(PLASMAC_HOME)/src/poisson_solver \
	$(PLASMAC_HOME)/src/rte_solver \
	$(PLASMAC_HOME)/src/sigma_solver \
	$(PLASMAC_HOME)/geometries/rod_sphere\
	$(PLASMAC_HOME)/physics/ito_plasma \
	$(PLASMAC_HOME)/physics/ito_plasma/cell_taggers/ito_plasma_streamer_tagger\
	$(PLASMAC_HOME)/physics/ito_plasma/time_steppers/ito_plasma_godunov\
	$(PLASMAC_HOME)/physics/ito_plasma/plasma_models/ito_plasma_air3


include $(CHOMBO_HOME)/mk/Make.example

