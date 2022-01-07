discharge-lib: source geometries 

physics: source advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer

chombo-lib:
	$(MAKE) --directory=$(CHOMBO_HOME) USE_EB=TRUE USE_MF=TRUE lib

all: discharge-lib physics

discharge-source: lib-chombo
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source 

discharge-geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries

advectiondiffusion: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion

brownianwalker: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker

cdrplasma: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma

electrostatics: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics

geometry: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Geometry

radiativetransfer: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/RadiativeTransfer

libclean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries pristine

allclean: libclean
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma          pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Geometry           pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/RadiativeTransfer  pristine
	$(RM) $(DISCHARGE_HOME)/Lib/*.a

pristine: allclean
	$(MAKE) --directory=$(CHOMBO_HOME) realclean

.PHONY: discharge-lib physics all source geometries advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer libclean allclean pristine
