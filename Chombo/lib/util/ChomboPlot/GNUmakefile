# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

CHOMBO_HOME = ../..
include $(CHOMBO_HOME)/mk/Make.defs

#this has to come after make.defs
AUTONAME = TRUE

ebase := chomboPlot

NEEDS_X11 := TRUE
XTRALIBFLAGS += $(X11LIBFLAGS)
XTRACPPFLAGS += $(X11INCFLAGS)

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTools ArrayView BoxTools

usage:
	@echo "This makefile supports targets:"
	@echo "     all util run $(ebase) usage "

# The `all-test', `run-only', `test-only' and $(ebase) targets are defined in `Make.rules'

ifeq ($(HAS_X11_MOTIF)-$(MPI),TRUE-FALSE)
  _all_target  = all-test
  _test_target = test-only
  _run_target  = run-only
else
  _all_target  = _no_target
  _test_target = _no_target
  _run_target  = _no_target
endif

all: $(_all_target)

util: $(_test_target)

run: $(_run_target)

_no_target:
	@echo 'HAS_X11_MOTIF is (' $(HAS_X11_MOTIF) ')'
	@echo 'MPI is (' $(MPI) ')'
	@echo "Warning: $(ebase) cannot be built without X11/Motif or with MPI."

.PHONY: usage all util run _no_target


include $(CHOMBO_HOME)/mk/Make.rules
