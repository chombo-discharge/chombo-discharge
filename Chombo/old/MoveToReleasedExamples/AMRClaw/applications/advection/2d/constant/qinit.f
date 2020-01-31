c     =====================================================
      subroutine qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      double precision qin, qout, xc, yc, win
      integer i,j


      common /comq/ qin, qout

      do i=1-mbc,mx+mbc
         xc = xlower + (i-0.5d0)*dx
         do j=1-mbc,my+mbc
            yc = ylower + (j - 0.5d0)*dy
            call cellave2(xc-dx/2, yc-dy/2,dx,dy,win)
            q(i,j,1) = win*qin + (1.d0-win)*qout
         enddo
      enddo

   50 continue

      return
      end


      double precision function fdisc(x,y)
      implicit none

      double precision x,y, xbar, ybar, radius, r2

      common /com_q0/ xbar, ybar, radius

c     Set up circle of radius 0.2 in center of grid
      r2 = (x-xbar)**2 + (y-ybar)**2
      if (r2 .le. radius**2) then
         fdisc = -1.d0
      else
         fdisc = 1.d0
      endif

      return

      end
