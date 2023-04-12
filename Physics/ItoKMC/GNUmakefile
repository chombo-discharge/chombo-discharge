# Contains various definitions of which Chombo libraries we need, and where the various chombo-discharge source files are located
include $(DISCHARGE_HOME)/Lib/Definitions.make

# This is for Chombo -- it will compile the files in src_dirs into a library. Chombo's makefiles will expand
# all files found in src_dirs into $objects so we can compile. 
src_dirs = $(ITOKMC_DIRS)
include $(CHOMBO_HOME)/mk/Make.example

# Name of the library
SrcLib = $(DISCHARGE_HOME)/Lib/lib$(ITOKMC_LIB)$(config).a

.DEFAULT_GOAL=ItoKMCPhysics

XTRACPPFLAGS += -I$(DISCHARGE_HOME)/Submodules/json/include

# For making the library.
ItoKMCPhysics: $(objects)
	if test -f $(SrcLib); then $(RM) $(SrcLib); fi
	ar -cr $(SrcLib) o/$(config)/*.o

# For cleaning up after ourselves. 
pristine:
	$(RM) -r [odpf]/$(config)
	if test -f $(SrcLib); then $(RM) $(SrcLib); fi

.PHONY: pristine ItoKMCPhysics
