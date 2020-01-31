      subroutine setprob
      implicit none

      double precision rhol, cl, rhor, cr
      double precision pi

      common /comaux/ rhol,cl,rhor,cr
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c
c     # Set the material parameters for the acoustic equations
c
      open(unit=7,file='setprob.data',status='old',form='formatted')

c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right:

      read(7,*) rhol
      read(7,*) cl
      read(7,*) rhor
      read(7,*) cr

      return
      end
