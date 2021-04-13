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
	$(PLASMAC_HOME)/src/cell_tagger \
	$(PLASMAC_HOME)/src/driver \
	$(PLASMAC_HOME)/src/geometry \
	$(PLASMAC_HOME)/src/global \
	$(PLASMAC_HOME)/src/particle \
	$(PLASMAC_HOME)/physics/geometry \
	$(PLASMAC_HOME)/geometries/tesselation\


include $(CHOMBO_HOME)/mk/Make.example

