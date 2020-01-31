c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # radially-symmetric pressure disturbance
c
       implicit none

       integer maxmx,maxmy, meqn, mbc, mx,my, maux
       double precision xlower, ylower, dx, dy
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision xi, yj, r, p0
       integer i, j
c
       do i = 1,mx
          xi = xlower + (i-0.5d0)*dx
          do j = 1,my
             yj = ylower + (j-0.5d0)*dy
             r = dsqrt((xi-0.25d0)**2 + (yj-0.4d0)**2)
             q(i,j,1) = p0(r)
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
          enddo
       enddo
       return
       end
