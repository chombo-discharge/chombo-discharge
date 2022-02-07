all: discharge-source discharge-geometries discharge-physics

chombo:
	$(MAKE) --directory=$(CHOMBO_HOME) USE_EB=TRUE USE_MF=TRUE lib

discharge-lib: discharge-source discharge-geometries discharge-physics

discharge-source: chombo
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source

discharge-geometries: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries

discharge-physics: discharge-source advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer

advectiondiffusion: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion

brownianwalker: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker

cdrplasma: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma

electrostatics: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics

geometry: discharge-source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Geometry

radiativetransfer: discharge-source
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
	find . -type f -name "*.ex" -delete
	find . -type f -name "*.d" -delete
	find . -type f -name "*.o" -delete

.PHONY: all chombo discharge-lib discharge-physics discharge-source discharge-geometries advectiondiffusion brownianwalker cdrplasma electrostatics geometry radiativetransfer libclean allclean pristine
