      program driver
c
c  Generic driver routine for claw2
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c

c     # This file is the driver routine for a clawpack example.
c     # To create the executable using Clawpack libraries files
c     # instead of Chombo libraries files, see Makefile

      implicit none


      integer maxmx, maxmy, mwork, mbc, meqn, mwaves, maux
      parameter (maxmx =   400)
      parameter (maxmy =   400)
      parameter (mwork = 1000000)

      parameter (mbc = 2)
      parameter (meqn = 4)
      parameter (mwaves = 3)
      parameter (maux = 7)
c       # NOTE: if maux>0 you must declare aux properly below!

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      double precision work(mwork)
      integer  mthlim(mwaves)

      integer maxlevel, level, refratios(20)
      common /comlevel/ maxlevel, level, refratios

      maxlevel = 0
      level = 0
c
c
      call claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop
      end
