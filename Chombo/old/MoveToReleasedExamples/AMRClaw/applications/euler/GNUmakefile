# -*- Mode: Makefile -*- 

%:
	@echo Doing nothing for AMRClaw/applications/euler in $(DIM)d for target $@

#XXX### This makefile requires GNU make.
#XXX### It builds and run the AMRNodeElliptic examples.
#XXX### This file is essentially the same as Chombo/example/GNUmakefile.
#XXX
#XXX### User-definable variables are in UPPER case.
#XXX### Internal variables are in lower case.
#XXX###
#XXX### Usage:
#XXX###   target:    purpose:
#XXX###   -------    ----------
#XXX###   all        builds all the executables, including the libraries they depend on
#XXX###   test       builds all the executables, but not the libraries
#XXX###   run        runs all the executables
#XXX###   clean      deletes files for this configuration
#XXX###   realclean  deletes files for all configurations that weren't in the original distribution
#XXX###
#XXX### Method:
#XXX###   The way this makefile works is the target specified by the user
#XXX###   on the command line is treated as the `action' to perform on all
#XXX###   the subdirectories in $dirs.  In the rules, the subdirectories are
#XXX###   treated like targets and the action is specified by a variable.
#XXX###   The makefile builds each directory target by invoking make
#XXX###   recursively on the directory with the action as the recursive target.
#XXX###   This assumes that the makefiles in the subdirectories implement all
#XXX###   the actions in the $(_all_actions) variable.
#XXX###
#XXX### Notes:
#XXX###
#XXX##################################################################################
#XXX
#XXXexport makefiles += example_AMRClaw_euler
#XXX
#XXX## This variable specifies where we are.
#XXXexport CHOMBO_HOME = ../../../../../lib
#XXX
#XXXinclude $(CHOMBO_HOME)/mk/Make.defs.defaults
#XXX
#XXX## This variable defines the names of the subdirectories to be built by default.
#XXX## One or more can be used as a target on the command line.
#XXXifeq ($(DIM),2)
#XXX  dirs = $(DIM)d/cylinder $(DIM)d/rp
#XXXelse
#XXX  dirs = $(DIM)d/rp $(DIM)d/rphex $(DIM)d/shockbubble $(DIMd/simplewave
#XXXendif
#XXX
#XXX##
#XXX## These are the actions (targets for the recursive make) supported by this
#XXX## makefile.  It assumes the makefiles in the subdirectories support them too.
#XXX##
#XXX_all_actions = all example example-only run clean
#XXX
#XXX.PHONY: $(_all_actions) $(dirs) usage vars
#XXX
#XXX##
#XXX## Default target
#XXX##
#XXXusage:
#XXX	@echo "This makefile supports targets:"
#XXX	@echo "     $(_all_actions) $(dirs) realclean usage vars"
#XXX
#XXX##
#XXX## For each action, set the action variable to tell the recursive
#XXX## makefile what to do.  There should be a separate dependency for
#XXX## every value in $(_all_actions).
#XXX##[NOTE: the `NODEPENDS=TRUE' variable disables automatic
#XXX##       dependency analysis in the recursive makefile, which is only
#XXX##       needed for actions that compile or link code.]
#XXX##
#XXX_action = all
#XXXall          : _action = all
#XXXexample      : _action = all
#XXXexample-only : _action = example-only
#XXXrun          : _action = run       NODEPENDS=TRUE
#XXXclean        : _action = clean     NODEPENDS=TRUE
#XXX
#XXX# this dependency ensures the action is done for every directory.
#XXX$(_all_actions) : $(dirs)
#XXX
#XXX# this rule performs the action
#XXX$(dirs):
#XXX	$(MAKE) --no-print-directory --directory $@ $(_action)
#XXX
#XXXrealclean:
#XXX	-find . -type d \( -name d -o -name p -o -name f -o -name o \) -depth -exec $(RM) -r {} \;
#XXX	-find . \( -name \*.ex -o -name \*.a -o -name \*_F.H -o -name core \) -exec $(RM) {} \;
#XXX	-find . \( -name \*.hdf5 -o -name \*.h5 -o -name \*~ -o -name work.pc\* -o -name gmon.out -o -name pout.\[0-9-\]\* -o -name leak.out\* -o -name .\#\* \) -exec $(RM) {} \;
#XXX
#XXX
#XXX# this rule prints out various 'make' variables to the screen for debugging
#XXXvars:
#XXX	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables $@
#XXX
