      subroutine setprob
      implicit none

      double precision qin,qout, ubar, vbar, xbar, ybar, radius

      common /comrp/ ubar,vbar
      common /comq/ qin, qout
      common /com_q0/ xbar, ybar, radius

c     # Set the velocity for scalar advection
c     # These values are passed to the Riemann solvers rpn2.f and rpt2.f
c     # in a common block

      open(unit=7,file='setprob.data',status='old',form='formatted')

      read(7,*) qin
      read(7,*) qout
      read(7,*) ubar
      read(7,*) vbar
      read(7,*) xbar
      read(7,*) ybar
      read(7,*) radius
      close(7)

      return
      end
