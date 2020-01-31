c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y + w*q_z = 0
c     # where u and v are a given velocity field.
c
c       -----------------------------------------------------------
c     # In advective form, with interface velocities specified in
c     # the auxiliary variable
c     # auxl(i,ma) contains auxiliary data for cells along this slice,
c     #   ma=1:   u-velocity at left edge of cell
c     #   ma=2:   v-velocity at bottom edge of cell
c     #   ma=3:   w-velocity at back edge of cell
c       -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  Note that in this advective form, the sum of
c     # amdq and apdq is not equal to a difference of fluxes except in the
c     # case of constant velocities.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
      implicit none

      integer ixyz, maxm, meqn, mwaves, mbc, mx, maux
c
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, maux)
      double precision auxr(1-mbc:maxm+mbc, maux)

      integer i


c     # Set wave, speed, and flux differences:
c     ------------------------------------------

      do i = 2-mbc, mx+mbc
         wave(i,1,1) = ql(i,1) - qr(i-1,1)
         s(i,1) = auxl(i,ixyz)
c        # The flux difference df = s*wave all goes in the downwind direction:
         amdq(i,1) = dmin1(auxl(i,ixyz), 0.d0) * wave(i,1,1)
         apdq(i,1) = dmax1(auxl(i,ixyz), 0.d0) * wave(i,1,1)
      enddo

      return
      end
