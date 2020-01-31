c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx,my,mz,maux
      double precision xlower, ylower, zlower, dx,dy,dz
c
      double precision  q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, meqn)
      double precision  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, maux)

      integer i,j,k
      double precision xi


       do k = 1,mz
          do j = 1,my
             do i = 1,mx
                xi = xlower + (i-0.5)*dx
                if (xi .lt. 0.5d0) then
                   q(i,j,k,1) = 1.d0
                else
                   q(i,j,k,1) = 0.d0
                endif
             enddo
          enddo
      enddo

      return
      end
