c     =====================================================
       subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                   xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit none

       integer maxmx, maxmy, maxmz,meqn,mbc,mx, my, mz, maux
       double precision xlower,ylower, zlower, dx, dy, dz

       double precision  q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &       1-mbc:maxmz+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &       1-mbc:maxmz+mbc, maux)


       double precision x0, y0, z0, xi, yj, zk, r, p0
       integer i,j,k

       x0 = 0.25d0
       y0 = 0.25d0
       z0 = 0.25d0

       do i = 1,mx
          xi = xlower + (i-0.5d0)*dx
          do j = 1,my
             yj = ylower + (j-0.5d0)*dy
             do k = 1,mz
                zk = zlower + (k-0.5d0)*dz
                r = sqrt((xi - x0)**2 + (yj - y0)**2 +
     &                (zk - z0)**2)
                q(i,j,k,1) = p0(r)
                q(i,j,k,2) = 0.d0
                q(i,j,k,3) = 0.d0
                q(i,j,k,4) = 0.d0
                enddo
             enddo
          enddo

      return
      end
