      subroutine setprob
      implicit none

      double precision z1, c1, z2, c2, pi

      common /compi/ pi
      common /comaux/ z1,c1,z2,c2

      pi = 4.d0*atan(1.d0)

      open(7,file='setprob.data')

c     # Piecewise constant medium with single interface as specified
c     # in setaux.f

c     # Impedance and sound speed in the two materials:

      read(7,*) z1
      read(7,*) c1
      read(7,*) z2
      read(7,*) c2
      close(7)

      return
      end
