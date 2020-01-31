# -*- Mode: Makefile -*- 

override DIM = 2

####define the location of Chombo
CHOMBO_HOME = ../../../../../../../lib
CCLAW_HOME  = ../../../..

# -----------------------------------------------------------------
# Include 2d f77 source files from current directory
F_SOURCES = \
          qinit.f \
          setprob.f \
	  estimate_error2.f \

# Include 2d f77 Riemann solvers for this directory
VPATH_LOCAL += ../rp
F_SOURCES += rpn2.f rpt2.f

# Include 2d f77 source files from CCLAW_HOME/src2d
VPATH_LOCAL += $(CCLAW_HOME)/src2d
F_SOURCES += \
             ClawPatch2.f \
             src2.f \
	     setaux.f \
             b4step2.f \
             step2.f \
             flux2.f  \
             qad2.f \
             bc2.f \
             fixcapaq2.f \
             average2.f \
             cellave2.f


# Include 2d/3d C++ source files from CCLAW_HOME/src
VPATH_LOCAL += $(CCLAW_HOME)/src
INCLUDES_LOCAL += -I$(CCLAW_HOME)/src -I$(CCLAW_HOME)/src2d

#[NOTE: Make.rules automatically gets $(ebase).cpp
CXX_SOURCES += AMRLevelClaw.cpp           \
               AMRLevelClawFactory.cpp    \
               ClawPatch.cpp              \
               CoarseAverageClaw.cpp      \
               FineInterpClaw.cpp

F_SOURCES += limiter.f philim.f

# gfortran needs to be told how to find 'claw.i'
override XTRAFFLAGS += $(INCLUDES_LOCAL)

# -----------------------------------------------------------------------
##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTimeDependent AMRTools BoxTools

ebase := amrClaw

INPUT = adv.inputs

include $(CHOMBO_HOME)/mk/Make.example.old
