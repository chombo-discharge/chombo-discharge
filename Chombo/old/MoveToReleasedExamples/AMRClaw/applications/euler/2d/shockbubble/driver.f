      program driver
c
c  Generic driver routine for claw2
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c
      implicit none

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c

      integer maxmx, maxmy, mwork, mbc, meqn, mwaves, maux

      parameter (maxmx =   640)
      parameter (maxmy =   160)
      parameter (mwork = 1271420)

      parameter (mbc = 2)
      parameter (meqn = 5)
      parameter (mwaves = 5)
      parameter (maux = 1)

c       # NOTE: if maux>0 you must declare aux properly below!

      double precision    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

      double precision  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision mthlim(mwaves)
      double precision  work(mwork)

c
      call claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop
      end
