# -*- Mode: Makefile -*- 

DIM = 2

#define the location of Chombo
CHOMBO_HOME = ../../../../../../../lib
CCLAW_HOME  = ../../../..


# ----------------------------------------------------------------------

# Include 2d f77 source files from current directory
F_SOURCES = \
	  qinit.f \
	  setprob.f \
	  estimate_error2.f

# Include 2d f77 Riemann solvers for this directory
VPATH_LOCAL += ../rp
F_SOURCES += rpn2ac.f rpt2ac.f

# Include 2d f77 source files from CCLAW_HOME/src2d
VPATH_LOCAL += $(CCLAW_HOME)/src2d
F_SOURCES += \
             ClawPatch2.f \
             src2.f \
             step2.f \
	     bc2.f \
             flux2.f  \
             b4step2.f \
             qad2.f \
	     setaux.f \
             average2.f \
             fixcapaq2.f


# Include 2d/3d F77 source files from CCLAW_HOME/src
VPATH_LOCAL += $(CCLAW_HOME)/src
INCLUDES_LOCAL += -I$(CCLAW_HOME)/src -I$(CCLAW_HOME)/src2d
F_SOURCES += limiter.f philim.f

#[NOTE: Make.rules automatically gets $(ebase).cpp
CXX_SOURCES += AMRLevelClaw.cpp           \
               AMRLevelClawFactory.cpp    \
               ClawPatch.cpp              \
               CoarseAverageClaw.cpp      \
               FineInterpClaw.cpp

# gfortran needs to be told how to find 'claw.i'
override XTRAFFLAGS += $(INCLUDES_LOCAL)

# -----------------------------------------------------------------------
##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTimeDependent AMRTools BoxTools

ebase := amrClaw

INPUT = wave.inputs

include $(CHOMBO_HOME)/mk/Make.example.old
