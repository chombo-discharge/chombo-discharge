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

      integer maxmx, maxmy, maxmz, mbc, mx, my, mz,maux
      double precision xlower, ylower, zlower, dx, dy, dz

      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)

      integer i, j, k
      double precision ubar, vbar, wbar

      common /comvel/ ubar, vbar, wbar

      do  k = 1-mbc,mz+mbc
	 do j = 1-mbc,my+mbc
	    do i = 1-mbc,mx+mbc
	       aux(i,j,k,1) = ubar
	       aux(i,j,k,2) = vbar
	       aux(i,j,k,3) = wbar
            enddo
         enddo
      enddo

      return
      end
