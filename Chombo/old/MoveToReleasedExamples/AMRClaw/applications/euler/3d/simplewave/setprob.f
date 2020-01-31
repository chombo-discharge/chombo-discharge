      subroutine setprob
      implicit none
      double precision gamma, gamma1, pi

      common /param/ gamma,gamma1
      common /compi/ pi

      open(unit=7,file='setprob.data',status='old',form='formatted')

      read(7,*) gamma
      gamma1 = gamma-1.d0
      close(7)

      pi = 4.d0*atan(1.d0)


      return
      end
