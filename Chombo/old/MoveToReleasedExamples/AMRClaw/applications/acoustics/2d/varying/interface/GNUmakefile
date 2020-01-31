# -*- Mode: Makefile -*- 

override DIM = 2

####define the location of Chombo
CHOMBO_HOME = ../../../../../../../../lib
CCLAW_HOME  = ../../../../..

# ----------------------------------------------------------------------

# Include 2d f77 source files from current directory
F_SOURCES = \
          setaux.f \
          qinit.f \
          setprob.f \
          p0.f

# Include 2d f77 Riemann solvers for this directory
VPATH_LOCAL += ../../rp
F_SOURCES += rpn2acv.f rpt2acv.f

# Include 2d f77 source files from CCLAW_HOME/src2d
VPATH_LOCAL += $(CCLAW_HOME)/src2d
F_SOURCES += \
             ClawPatch2.f \
             src2.f \
             b4step2.f \
             step2.f \
             flux2.f  \
             qad2.f \
             bc2.f \
             fixcapaq2.f \
             average2.f \
             estimate_error2.f

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


#ifeq ($(fname),gfortran)
# -I$(CCLAW_HOME)/src2d is used here as a temporary fix -- 
#  appears gfortran (new with gcc4.0) needs to be told explicitly where to find include files.
#  -ffixed-line-length-none is to avoid a file truncation warning (only happened in the comments)
#  but i don't know how to turn on flags for a certain compiler...
override XTRAFFLAGS += $(INCLUDES_LOCAL) #-ffixed-line-length-none
#endif

# -----------------------------------------------------------------------
##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTimeDependent AMRTools BoxTools

ebase := amrClaw

INPUT = acv.inputs.TEST

include $(CHOMBO_HOME)/mk/Make.example.old
