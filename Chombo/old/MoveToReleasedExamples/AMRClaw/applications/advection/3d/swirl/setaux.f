c     ==================================================================
      subroutine setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c
c     # advection
c     #    aux(i,j,k,1) is u velocity on left face of cell
c     #    aux(i,j,k,2) is v velocity on bottom face of cell
c     #    aux(i,j,k,3) is w velocity on back face of cell
c
c
      implicit none

      integer maxmx, maxmy, maxmz, mx, my, mz, mbc, maux, i,j,k,m
      double precision xlower, ylower, zlower, dx, dy, dz,
     &      pi, pi2, dx2, dy2, dz2, compute_u, compute_v, compute_w
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)

      double precision xc, yc, zc

      dx2 = 0.5d0*dx
      dy2 = 0.5d0*dy
      dz2 = 0.5d0*dz

      do  k = 1-mbc,mz+mbc
         zc = zlower + (k-0.5d0)*dz
	 do j = 1-mbc,my+mbc
            yc = ylower + (j-0.5d0)*dy
	    do i = 1-mbc,mx+mbc
               xc = xlower + (i-0.5)*dx
	       aux(i,j,k,1) = compute_u(xc -dx2, yc ,zc)
	       aux(i,j,k,2) = compute_v(xc,yc-dy2, zc)
	       aux(i,j,k,3) = compute_w(xc,yc, zc-dz2)
            enddo
         enddo
      enddo

      return
      end

      double precision function compute_u(x,y,z)
      implicit none
      double precision x,y,z,pi, pi2

      common /compsi/ pi

      pi2 = 2.d0*pi
      compute_u = 2.d0 * dsin(pi*x)**2 * dsin(pi2*y) * dsin(pi2*z)
      end

      double precision function compute_v(x,y,z)
      implicit none
      double precision x,y,z,pi, pi2

      common /compsi/ pi

      pi2 = 2.d0*pi
      compute_v = -dsin(pi2*x) * dsin(pi*y)**2 * dsin(pi2*z)
      end

      double precision function compute_w(x,y,z)
      implicit none
      double precision x,y,z,pi, pi2

      common /compsi/ pi

      pi2 = 2.d0*pi
      compute_w = -dsin(pi2*x) * dsin(pi2*y) * dsin(pi*z)**2
      end
