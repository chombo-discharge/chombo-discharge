      subroutine setprob
      implicit none

      double precision rho, bulk, cc, zz, width, pi

      common /cparam/ rho,bulk,cc,zz
      common /initq/ width
      common /compi/ pi
c
c     # Set the material parameters for the acoustic equations
c     # Passed to the Riemann solver rp1.f in a common block
c
      open(7, file='setprob.data')

c     # density:
      read(7,*) rho

c     # bulk modulus:
      read(7,*) bulk

      read(7,*) width

c     # sound speed:
      cc = dsqrt(bulk/rho)

c     # impedance:
      zz = cc*rho

      close(7)

      return
      end
