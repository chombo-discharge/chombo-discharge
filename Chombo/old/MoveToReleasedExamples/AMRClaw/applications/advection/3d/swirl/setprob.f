      subroutine setprob
      implicit none

      double precision pi, tperiod, pi2

      common /compi/ pi
      common /comvt/ tperiod,pi2

      open(unit=7,file='setprob.data',status='old',form='formatted')
      read(7,*) tperiod
      close(7)

c     # compute pi, used in psi.f
      pi = 4.d0 * datan(1.d0)

c     # save 2*pi and tperiod in common block for use in b4step2:
      pi2 = 2.d0*pi


      return
      end
