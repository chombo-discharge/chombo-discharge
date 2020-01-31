c
c
c =========================================================
       subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my, mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # smooth density perturbation of magnitude rhoamp centered at x0
c     # with p, u chosen so it is a simple 3-wave before shock formation
c
c
      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx, my, mz, maux

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc,maux)

      double precision xlower, ylower, zlower, dx, dy, dz
      double precision gamma, gamma1

      integer k,j,i
      double precision rhoamp, x0, xcell, ycell, zcell, xp, yp, zp
      double precision pres, u, v, w, velnorm
      common /param/ gamma, gamma1
c
c
      rhoamp = 0.25d0
      x0 = -1.75d0
      do k = 1-mbc,mz+mbc
         zcell = zlower + (k-0.5d0)*dz
         do j = 1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
            do i = 1-mbc,mx+mbc
               xcell = xlower + (i-0.5d0)*dx
               call mapc2p(xcell,ycell,zcell,xp,yp,zp);
               q(i,j,k,1) = 1.d0 + rhoamp*dexp(-5.d0*(xp-x0)**2)
               pres = q(i,j,k,1)**gamma
               u = 2.d0/gamma1 * (dsqrt(gamma*q(i,j,k,1)**gamma1) -
     &            dsqrt(gamma)) + 2.d0
               v = 0.d0
               w = 0.d0
               q(i,j,k,2) = q(i,j,k,1)*u
               q(i,j,k,3) = q(i,j,k,1)*v
               q(i,j,k,4) = q(i,j,k,1)*w
               velnorm = u*u + v*v + w*w
               q(i,j,k,5) = pres / gamma1 + 0.5d0*q(i,j,k,1)*velnorm
            enddo
         enddo
      enddo
c
         write(6,*) 'Done with qinit'

      return
      end
