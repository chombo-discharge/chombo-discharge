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
	$(DISCHARGE_HOME)/Source/ImplicitFunctions \
	$(DISCHARGE_HOME)/Source/Multifluid \
	$(DISCHARGE_HOME)/Source/Utilities \
	$(DISCHARGE_HOME)/Source/Driver \
	$(DISCHARGE_HOME)/Source/Particle \
	$(DISCHARGE_HOME)/Source/ItoDiffusion \
	$(DISCHARGE_HOME)//Source/Electrostatics \
	$(DISCHARGE_HOME)/Source/RadiativeTransfer \
	$(DISCHARGE_HOME)/Source/SigmaSolver \
	$(DISCHARGE_HOME)/Geometries/RodDielectric\
	$(DISCHARGE_HOME)/Physics/ItoPlasma \
	$(DISCHARGE_HOME)/Physics/ItoPlasma/CellTaggers/ito_plasma_streamer_tagger\
	$(DISCHARGE_HOME)/Physics/ItoPlasma/TimeSteppers/ito_plasma_godunov\
	$(DISCHARGE_HOME)/Physics/ItoPlasma/PlasmaModels/ito_plasma_air3


include $(CHOMBO_HOME)/mk/Make.example

