include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system. 
ebase    = program
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebaseobject): dependencies
.DEFAULT_GOAL=$(ebase)

# Build dependencies if they do not exis. 
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME) discharge-lib
	$(MAKE) --directory=$(DISCHARGE_HOME) cdrplasma

# Make advection-diffusion headers and library visible. 
XTRACPPFLAGS += $(CDRPLASMA_INCLUDE) -I$(DISCHARGE_HOME)/Submodules/json/include
XTRALIBFLAGS += $(addprefix -l, $(CDRPLASMA_LIB))$(config)
