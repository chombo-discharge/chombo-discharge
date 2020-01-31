#
# This version of the standard Clawpack Makefile is the version used with ChomboClaw.
# The purpose is to give users who wish to migrate their code to ChomboClaw
# an easy way to run a single grid calculation using Clawpack.
#
# NOTE : The Clawpack source code should be installed and built and available under
# the directory indicated by the environment variable CLAW.
#
#               Makefile for the clawpack code:
#
#       To make an executable, type:            gmake -f Makefile
#       (from the application sub-directory)
#
#
#  This file is to be used with ChomboClaw for testing purposes.
#
# Donna Calhoun; 4/6/2006.
#


FFLAGS = -O
LFLAGS = $(FFLAGS)
F77    = g77 -c
LINK   = g77

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  driver.o \
  qinit.o \
  ../rp/rpn2eu5.o \
  ../rp/rpt2eu5.o \
  setprob.o \
  src2.o \
  setaux.o \
  ../../../../src2d/cellave2.o \
  fdisc.o

LIBOBJECTS = \
  $(CLAW)/clawpack/2d/lib/claw2ez.o \
  $(CLAW)/clawpack/2d/lib/claw2.o \
  $(CLAW)/clawpack/2d/lib/b4step2.o \
  $(CLAW)/clawpack/2d/lib/bc2.o \
  $(CLAW)/clawpack/2d/lib/step2.o \
  $(CLAW)/clawpack/2d/lib/step2ds.o \
  $(CLAW)/clawpack/2d/lib/dimsp2.o \
  $(CLAW)/clawpack/2d/lib/flux2.o \
  $(CLAW)/clawpack/2d/lib/copyq2.o \
  $(CLAW)/clawpack/2d/lib/limiter.o \
  $(CLAW)/clawpack/2d/lib/philim.o \
  $(CLAW)/clawpack/2d/lib/out2.o \


xclaw: $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclaw


clean:
	rm -f *.o core* xclaw

realclean:
	rm -f *.o core* xclaw fort.*

### DO NOT remove this line - make depends on it ###
