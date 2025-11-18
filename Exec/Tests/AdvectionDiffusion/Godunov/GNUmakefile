include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system. 
ebase    = main
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebaseobject): dependencies
.DEFAULT_GOAL=$(ebase)

# Build dependencies.
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME) discharge-lib
	$(MAKE) --directory=$(DISCHARGE_HOME) advectiondiffusion

# Make advection-diffusion headers and library visible. 
XTRACPPFLAGS += $(ADVDIFF_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(ADVDIFF_LIB))$(config)
