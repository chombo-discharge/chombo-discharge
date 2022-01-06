lib: source geometries 

physics: source advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer

all: lib physics

source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source 

geometries: source
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

pristine: allclean
	$(RM) $(DISCHARGE_HOME)/Lib/*.a


.PHONY: lib physics all source geometries advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer libclean allclean pristine
