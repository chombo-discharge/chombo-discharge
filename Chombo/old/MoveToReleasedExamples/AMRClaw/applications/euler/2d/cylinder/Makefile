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
#  This file is to be used with ChomboClaw for testing purposes.
#
# Donna Calhoun; 4/6/2006.
#

CHOMBO_HOME = ../../../../src2d

FFLAGS = -O
LFLAGS = $(FFLAGS)
F77    = ${FC} -c
LINK   = ${LD}

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  driver.o \
  qinit.o \
  ../rp/rpn2euq.o \
  ../rp/rpt2euq.o \
  setprob.o \
  bc2-curve.o \
  mapc2p.o \
  setaux.o \
  ${CHOMBO_HOME}/setquadinfo.o \
  ${CHOMBO_HOME}/getquadinfo.o \
  ${CHOMBO_HOME}/cellave2.o

LIBOBJECTS = \
  $(CLAW)/clawpack/2d/lib/claw2ez.o \
  $(CLAW)/clawpack/2d/lib/claw2.o \
  $(CLAW)/clawpack/2d/lib/src2.o \
  $(CLAW)/clawpack/2d/lib/b4step2.o \
  $(CLAW)/clawpack/2d/lib/step2.o \
  $(CLAW)/clawpack/2d/lib/step2ds.o \
  $(CLAW)/clawpack/2d/lib/dimsp2.o \
  $(CLAW)/clawpack/2d/lib/flux2.o \
  $(CLAW)/clawpack/2d/lib/copyq2.o \
  $(CLAW)/clawpack/2d/lib/limiter.o \
  $(CLAW)/clawpack/2d/lib/philim.o \
  $(CLAW)/clawpack/2d/lib/out2.o

xclaw: $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) -o xclaw

clean :
	rm -f *.o core.* xclaw ../rp/*.o

realclean :
	rm -f *.o core.* xclaw ../rp/*.o fort.*

### DO NOT remove this line - make depends on it ###
