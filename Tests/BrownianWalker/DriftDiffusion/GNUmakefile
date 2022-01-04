include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system. 
ebase    = main
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebase): dependencies

# Build dependencies if they do not exis. 
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME) lib
	$(MAKE) --directory=$(DISCHARGE_HOME) brownianwalker

# Make advection-diffusion headers and library visible. 
XTRACPPFLAGS += $(BROWNIAN_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(BROWNIAN_LIB))$(config)
