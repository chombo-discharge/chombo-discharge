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
	$(DISCHARGE_HOME)/Source/CellTagger \
	$(DISCHARGE_HOME)/Source/Driver \
	$(DISCHARGE_HOME)/Source/Elliptic \
	$(DISCHARGE_HOME)/Source/geometry \
	$(DISCHARGE_HOME)/Source/rte_solver \
	$(DISCHARGE_HOME)/Source/sigma_solver \
	$(DISCHARGE_HOME)/Source/global \
	$(DISCHARGE_HOME)/Source/FieldSolver \
	$(DISCHARGE_HOME)/Source/particle \
	$(DISCHARGE_HOME)/geometries/rod_plane_profile\
	$(DISCHARGE_HOME)/Physics/Electrostatics \


include $(CHOMBO_HOME)/mk/Make.example

