source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source lib

geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries lib

advectiondiffusion: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion lib

brownianwalker: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker lib

cdrplasma: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma lib

electrostatics: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics lib

geometry: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Geometry lib

radiativetransfer: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/RadiativeTransfer lib

clean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries pristine

realclean: clean
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma          pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Geometry           pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/RadiativeTransfer  pristine

physics: source advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer

lib: source geometries 

all: lib physics

