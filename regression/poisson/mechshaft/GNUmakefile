include $(CHOMBO_HOME)/mk/Make.defs

USE_EB=TRUE
USE_MF=TRUE
DIM=3


# Base file containing int main()
ebase := main

LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools


# Target
all: all-test

base_dir = .
src_dirs = $(DISCHARGE_HOME)/src \
	$(DISCHARGE_HOME)/src/amr_mesh \
	$(DISCHARGE_HOME)/src/cell_tagger \
	$(DISCHARGE_HOME)/src/driver \
	$(DISCHARGE_HOME)/src/elliptic \
	$(DISCHARGE_HOME)/src/geometry \
	$(DISCHARGE_HOME)/src/rte_solver \
	$(DISCHARGE_HOME)/src/sigma_solver \
	$(DISCHARGE_HOME)/src/global \
	$(DISCHARGE_HOME)/src/field_solver \
	$(DISCHARGE_HOME)/src/particle \
	$(DISCHARGE_HOME)/geometries/mechanical_shaft\
	$(DISCHARGE_HOME)/physics/poisson \


include $(CHOMBO_HOME)/mk/Make.example

