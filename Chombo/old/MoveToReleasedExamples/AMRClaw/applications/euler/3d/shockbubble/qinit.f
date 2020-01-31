c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx, my, mz,maux
      double precision xlower,ylower, zlower,dx,dy,dz

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,maux)

      double precision gamma, gamma1, qin(5),qout(5)
      double precision rinf, vinf, einf
      double precision x0,y0,z0,r0
      double precision xe,ye,ze,d2,win
      integer i,j,k,m

      common /param/ gamma, gamma1
      common /comic/ qin,qout
      common /cominf/ rinf,vinf,einf
      common /cdisc/ x0,y0,z0,r0

      do i = 1,mx
         xe = xlower + (i-1)*dx
         do j = 1,my
            ye = ylower + (j-1)*dy
            do k = 1,mz
               ze = zlower + (k - 1)*dz
               d2 = (xe - x0)**2 + (ye - y0)**2 +
     &               (ze - z0)**2
               if (d2 .le. (r0+2*dx)**2) then
                  call cellave3(xe,ye,ze,dx,dy,dz,win)
               else
                  win = 0.0
               endif
               do m=1,meqn
                  q(i,j,k,m) = win*qin(m) + (1.d0-win)*qout(m)
               enddo

               if (xe < 0.2d0) then
c                 # behind shock:
                  q(i,j,k,1) = rinf
                  q(i,j,k,2) = rinf*vinf
                  q(i,j,k,3) = 0.d0
                  q(i,j,k,4) = 0.d0
                  q(i,j,k,5) = einf
               end if
            enddo
         enddo
      enddo

      return
      end

c     =================================================
      double precision function fdisc(x,y,z)
c     =================================================
      implicit none

      double precision x,y,z,x0,y0,z0,r0
      common/cdisc/ x0,y0, z0, r0

c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the
c     # left of the curve and positive to the right
c     # idisc specifies the nature of the discontinuity for two
c     # particular cases (a straight line and circle) but this routine
c     # can be modified for any other curve.
c

c     # circle of radius r0:

c     # fdisc < 0 inside sphere, and fdisc >= 0 outside the sphere

      fdisc = ((x-x0)**2 + (y-y0)**2 + (z-z0)**2)-r0**2
c
      return
      end
