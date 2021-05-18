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
src_dirs = $(DISCHARGE_HOME)/Source \
	$(DISCHARGE_HOME)/Source/AmrMesh \
	$(DISCHARGE_HOME)/Source/CdrSolver \
	$(DISCHARGE_HOME)/Source/CellTagger \
	$(DISCHARGE_HOME)/Source/Elliptic \
	$(DISCHARGE_HOME)/Source/geometry \
	$(DISCHARGE_HOME)/Source/global \
	$(DISCHARGE_HOME)/Source/Driver \
	$(DISCHARGE_HOME)/Source/particle \
	$(DISCHARGE_HOME)/Source/ito_solver \
	$(DISCHARGE_HOME)/Source/FieldSolver \
	$(DISCHARGE_HOME)/Source/rte_solver \
	$(DISCHARGE_HOME)/Source/sigma_solver \
	$(DISCHARGE_HOME)/Geometries/rod_dielectric\
	$(DISCHARGE_HOME)/Physics/ito_plasma \
	$(DISCHARGE_HOME)/Physics/ito_plasma/cell_taggers/ito_plasma_streamer_tagger\
	$(DISCHARGE_HOME)/Physics/ito_plasma/time_steppers/ito_plasma_godunov\
	$(DISCHARGE_HOME)/Physics/ito_plasma/plasma_models/ito_plasma_air3


include $(CHOMBO_HOME)/mk/Make.example

