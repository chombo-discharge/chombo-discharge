include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system.
ebase    = main
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code.
$(ebaseobject): dependencies
.DEFAULT_GOAL=$(ebase)

# Build dependencies if they do not exist.
dependencies:
	$(MAKE) --directory=$(DISCHARGE_HOME) discharge-lib
	$(MAKE) --directory=$(DISCHARGE_HOME) brownianwalker

# Make BrownianWalker headers/library visible, and the Dev SoA staging headers (validation only).
XTRACPPFLAGS += $(BROWNIAN_INCLUDE)
XTRACPPFLAGS += -I$(DISCHARGE_HOME)/Dev
XTRALIBFLAGS += $(addprefix -l, $(BROWNIAN_LIB))$(config)
