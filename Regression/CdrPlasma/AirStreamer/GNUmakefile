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
	$(DISCHARGE_HOME)//Source/ConvectionDiffusionReaction \
	$(DISCHARGE_HOME)/Source/CellTagger \
	$(DISCHARGE_HOME)/Source/Elliptic \
	$(DISCHARGE_HOME)/Source/Geometry \
	$(DISCHARGE_HOME)/Source/global \
	$(DISCHARGE_HOME)/Source/Driver \
	$(DISCHARGE_HOME)/Source/Particle \
	$(DISCHARGE_HOME)//Source/Electrostatics \
	$(DISCHARGE_HOME)/Source/RadiativeTransfer \
	$(DISCHARGE_HOME)/Source/SigmaSolver \
	$(DISCHARGE_HOME)/Geometries/rod_dielectric\
	$(DISCHARGE_HOME)/Physics/cdr_plasma \
	$(DISCHARGE_HOME)/Physics/cdr_plasma/cell_taggers/streamer_tagger\
	$(DISCHARGE_HOME)/Physics/cdr_plasma/time_steppers/imex_sdc\
	$(DISCHARGE_HOME)/Physics/cdr_plasma/plasma_models/air3_bourdon


include $(CHOMBO_HOME)/mk/Make.example

