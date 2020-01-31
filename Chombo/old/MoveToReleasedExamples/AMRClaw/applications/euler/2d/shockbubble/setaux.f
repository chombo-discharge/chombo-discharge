c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c     # aux(i,j,1) = y coordinate of cell center for cylindrical source terms
c
c
      implicit none

      integer maxmx, maxmy, mbc, mx,my, maux
      double precision xlower, ylower, dx, dy

      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i,j
      double precision rj

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            rj = ylower + (j-0.5d0)*dy
            aux(i,j,1) = rj
         enddo
      enddo

       return
       end
