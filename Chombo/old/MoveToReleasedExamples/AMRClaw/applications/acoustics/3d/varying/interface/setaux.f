c     ==================================================================
      subroutine setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c
c     # acoustics in a heterogeneous medium:
c     #  aux(i,j,k,1) = impedance Z in (i,j) cell
c     #  aux(i,j,k,2) = sound speed c in (i,j) cell
c
c
      implicit none

      integer maxmx, maxmy, maxmz, mbc, mx, my, mz, maux
      double precision xlower,ylower, zlower, dx, dy, dz
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)


      double precision z1, c1, z2, c2, xcell, ycell, zcell
      integer i, j, k

      common /comaux/ z1,c1,z2,c2


      do  k = 1-mbc,mz+mbc
         zcell = zlower + (k-0.5d0)*dz
	 do j = 1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
	    do i = 1-mbc,mx+mbc
               xcell = xlower + (i-0.5)*dx
               if (xcell .lt. 0.5d0) then
	          aux(i,j,k,1) = z1
	          aux(i,j,k,2) = c1
                else
	          aux(i,j,k,1) = z2
	          aux(i,j,k,2) = c2
                endif
            enddo
         enddo
      enddo

      return
      end
