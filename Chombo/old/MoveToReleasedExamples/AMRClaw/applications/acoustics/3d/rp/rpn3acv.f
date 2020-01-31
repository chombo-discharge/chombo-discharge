c
c
c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Riemann solver for the acoustics equations in 3d, with varying
c     # material properties.
c
c     # auxl(i,1) holds impedance Z,
c     # auxl(i,2) holds sound speed c,
c
c     # Note that although there are 4 eigenvectors, two eigenvalues are
c     # always zero and so we only need to compute 2 waves.
c
c     # Solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit none
      integer ixyz, maxm, meqn,mwaves, mbc, mx, maux
      double precision  wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision     s(1-mbc:maxm+mbc, mwaves)
      double precision    ql(1-mbc:maxm+mbc, meqn)
      double precision    qr(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc, maux)
      double precision  auxr(1-mbc:maxm+mbc, maux)


      integer mu,mv, mw, i, m
      double precision  delta(3), zi, zim, a1, a2, area

c     # set mu to point to  the component of the system that corresponds
c     # to velocity in the direction of this slice, mv to the orthogonal
c     # velocity.
c
c
      if (ixyz .eq. 1) then
          mu = 2
          mv = 3
          mw = 4
        else if (ixyz .eq. 2) then
          mu = 3
          mv = 4
          mw = 2
        else if (ixyz .eq. 3) then
          mu = 4
          mv = 2
          mw = 3
        endif

c
c     # split the jump in q at each interface into waves
c     # The jump is split into a leftgoing wave traveling at speed -c
c     # relative to the material properties to the left of the interface,
c     # and a rightgoing wave traveling at speed +c
c     # relative to the material properties to the right of the interface,
c
c     # find a1 and a2, the coefficients of the 2 eigenvectors:
      do i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,mu) - qr(i-1,mu)
c        # impedances:
         zi = auxl(i,1)
         zim = auxl(i-1,1)

         a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
         a2 =  (delta(1) + zim*delta(2)) / (zim + zi)


c
c        # Compute the waves.

         wave(i,1,1) = -a1*zim
         wave(i,mu,1) = a1
         wave(i,mv,1) = 0.d0
         wave(i,mw,1) = 0.d0
         s(i,1) = -auxl(i-1,2)

         wave(i,1,2) = a2*zi
         wave(i,mu,2) = a2
         wave(i,mv,2) = 0.d0
         wave(i,mw,2) = 0.d0
         s(i,2) = auxl(i,2)
      enddo

c     # compute the leftgoing and rightgoing flux differences:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      area = 1.0
      do m = 1,meqn
         do i = 2-mbc, mx+mbc
            s(i,1) = area*s(i,1)
            s(i,2) = area*s(i,2)
            amdq(i,m) = s(i,1)*wave(i,m,1)
            apdq(i,m) = s(i,2)*wave(i,m,2)
         enddo
      enddo

      end
