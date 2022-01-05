include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system. 
ebase    = main
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebaseobject): dependencies

# Build dependencies if they do not exis. 
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME)
	$(MAKE) --directory=$(DISCHARGE_HOME) electrostatics

# Make advection-diffusion headers and library visible. 
XTRACPPFLAGS += $(ELECTROSTATICS_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(ELECTROSTATICS_LIB))$(config)
