include $(DISCHARGE_HOME)/Lib/Definitions.make

# Things for the Chombo makefile system. 
ebase    = main
include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebaseobject): dependencies
.DEFAULT_GOAL=$(ebase)

# Build dependencies if they do not exist. This application uses both the CdrPlasma and ItoKMC modules.
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME) discharge-lib
	$(MAKE) --directory=$(DISCHARGE_HOME) cdrplasma
	$(MAKE) --directory=$(DISCHARGE_HOME) ito-kmc

# Make the CdrPlasma and ItoKMC headers and libraries visible. 
XTRACPPFLAGS += $(CDRPLASMA_INCLUDE) $(ITOKMC_INCLUDE) -I$(DISCHARGE_HOME)/Submodules/json/include
XTRALIBFLAGS += $(addprefix -l, $(CDRPLASMA_LIB))$(config) $(addprefix -l, $(ITOKMC_LIB))$(config)
