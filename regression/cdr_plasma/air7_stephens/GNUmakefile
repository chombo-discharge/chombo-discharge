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
src_dirs = $(DISCHARGE_HOME)/src \
	$(DISCHARGE_HOME)/src/AmrMesh \
	$(DISCHARGE_HOME)/src/CdrSolver \
	$(DISCHARGE_HOME)/src/CellTagger \
	$(DISCHARGE_HOME)/src/elliptic \
	$(DISCHARGE_HOME)/src/geometry \
	$(DISCHARGE_HOME)/src/global \
	$(DISCHARGE_HOME)/src/Driver \
	$(DISCHARGE_HOME)/src/particle \
	$(DISCHARGE_HOME)/src/FieldSolver \
	$(DISCHARGE_HOME)/src/rte_solver \
	$(DISCHARGE_HOME)/src/sigma_solver \
	$(DISCHARGE_HOME)/geometries/rod_dielectric\
	$(DISCHARGE_HOME)/physics/cdr_plasma \
	$(DISCHARGE_HOME)/physics/cdr_plasma/cell_taggers/streamer_tagger\
	$(DISCHARGE_HOME)/physics/cdr_plasma/time_steppers/godunov\
	$(DISCHARGE_HOME)/physics/cdr_plasma/plasma_models/air7_stephens


include $(CHOMBO_HOME)/mk/Make.example

