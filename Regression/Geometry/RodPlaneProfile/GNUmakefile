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
	$(DISCHARGE_HOME)/Source/Geometry \
	$(DISCHARGE_HOME)/Source/ImplicitFunctions \
	$(DISCHARGE_HOME)/Source/Multifluid \
	$(DISCHARGE_HOME)/Source/Utilities \
	$(DISCHARGE_HOME)/Source/Particle \
	$(DISCHARGE_HOME)/Physics/Geometry \
	$(DISCHARGE_HOME)/Geometries/RodPlaneProfile\


include $(CHOMBO_HOME)/mk/Make.example

