# -*- Mode: Makefile -*- 

override DIM = 2

####define the location of Chombo
CHOMBO_HOME = ../../../../../../../lib
CCLAW_HOME  = ../../../..

#  # needed to determine compiler commands
#  include $(CHOMBO_HOME)/mk/Make.defs
#
#  cxxname := $(notdir $(firstword $(CXX)))
#  fname   := $(notdir $(firstword $(FC)))

# ----------------------------------------------------------------------
# Include 2d f77 source files from current directory
F_SOURCES = \
          setaux.f \
          qinit.f \
          setprob.f \
          mapc2p.f \
          bc2-curve.f \
          estimate_error2.f

# Include 2d f77 Riemann solvers for this directory
VPATH_LOCAL += ../rpquad
F_SOURCES += rpn2euq.f rpt2euq.f

# Include 2d f77 source files from CCLAW_HOME/src2d
VPATH_LOCAL += $(CCLAW_HOME)/src2d
F_SOURCES += \
             ClawPatch2.f \
             src2.f \
             step2.f \
             flux2.f  \
             b4step2.f \
             qad2.f \
             average2.f \
             setquadinfo.f \
             getquadinfo.f \
             cellave2.f \
             fixcapaq2.f

# ----------------------------------------------------------------------
# Include 2d/3d source files from CCLAW_HOME/src
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

INPUT = shock.inputs.TEST

include $(CHOMBO_HOME)/mk/Make.example.old
