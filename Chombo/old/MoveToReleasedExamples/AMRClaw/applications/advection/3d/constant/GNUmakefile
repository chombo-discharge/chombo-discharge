# -*- Mode: Makefile -*- 

override DIM = 3

####define the location of Chombo
CHOMBO_HOME = ../../../../../../../lib
CCLAW_HOME  = ../../../..

# ------------------------------------------------------------------
# Include 3d f77 source files from current directory
F_SOURCES = \
          setaux.f \
          qinit.f \
          setprob.f

# Include 3d f7 Riemann solvers for this directory
VPATH_LOCAL += ../rp
F_SOURCES += rpn3adv.f rpt3adv.f rptt3adv.f

# Include 3d f77 source files from CCLAW_HOME/src3d
VPATH_LOCAL += $(CCLAW_HOME)/src3d
F_SOURCES += \
             ClawPatch3.f \
             src3.f \
             estimate_error3.f \
             average3.f \
             fixcapaq3.f \
             step3.f \
             flux3.f  \
             qad3.f \
             b4step3.f \
             bc3.f


# Include 2d/3d F77 source files from CCLAW_HOME/src
VPATH_LOCAL += $(CCLAW_HOME)/src
INCLUDES_LOCAL += -I$(CCLAW_HOME)/src -I$(CCLAW_HOME)/src3d

F_SOURCES += limiter.f philim.f

#[NOTE: Make.rules automatically gets $(ebase).cpp
CXX_SOURCES += AMRLevelClaw.cpp           \
               AMRLevelClawFactory.cpp    \
               CoarseAverageClaw.cpp      \
               FineInterpClaw.cpp         \
               ClawPatch.cpp

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
