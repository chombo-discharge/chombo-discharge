c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================

      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx, my,mz, maux
      double precision xlower, ylower, zlower, dx, dy, dz
c
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, maux)

      integer i, j, k
      double precision qin, qout, xbar, ybar, zbar, radius
      double precision x,y,z,r

      common /comq/ qin, qout
      common /comgeo/ xbar, ybar, zbar, radius

c
c     # set concentration profile
c     ---------------------------
c
       do i = 1,mx
          x = xlower + (i-0.5d0)*dx
          do j = 1,my
             y = ylower + (j-0.5d0)*dy
             do k = 1,mz
                z = zlower + (k-0.5d0)*dz
                r = (x-xbar)**2 + (y-ybar)**2 + (z-zbar)**2
                if (r .le. radius**2) then
                   q(i,j,k,1) = qin
                else
                   q(i,j,k,1) = qout
                endif
             enddo
          enddo
       enddo


      return
      end
