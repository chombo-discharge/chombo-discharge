      subroutine setprob
      implicit none

      double precision ubar, vbar, wbar, qin, qout
      double precision xbar,ybar, zbar, radius

      common /comvel/ ubar,vbar, wbar
      common /comq/ qin, qout
      common /comgeo/ xbar, ybar, zbar, radius
c
c     # Set the velocity for scalar advection
c     # These values are set in aux array in setaux.f and
c     # used in Riemann solvers.
c
      open(unit=7,file='setprob.data',status='old',form='formatted')

      read(7,*) qin
      read(7,*) qout
      read(7,*) ubar
      read(7,*) vbar
      read(7,*) wbar
      read(7,*) xbar
      read(7,*) ybar
      read(7,*) zbar
      read(7,*) radius
      close(7)

      return
      end
