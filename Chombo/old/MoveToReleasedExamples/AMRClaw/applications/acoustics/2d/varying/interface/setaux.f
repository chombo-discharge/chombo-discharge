c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c     # variable coefficient acoustics
c     #  aux(i,j,1) = density rho in (i,j) cell
c     #  aux(i,j,2) = sound speed c in (i,j) cell
c
c     # Piecewise constant medium with single interface at x=0.5
c     # Density and sound speed to left and right are set in setprob.f

c
c
      implicit none

      integer maxmx, maxmy, mbc, mx, my, maux
      double precision dx, dy, xlower, ylower
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      double precision rhol, cl, rhor, cr
      double precision xcell, ycell
      integer i,j

      common /comaux/ rhol,cl,rhor,cr
c

      do  j=1-mbc,my+mbc
         do  i=1-mbc,mx+mbc
            xcell = xlower + (i - 0.5d0)*dx
            ycell = ylower + (j - 0.5d0)*dy
            if (xcell .lt. 0.5d0) then
               aux(i,j,1) = rhol
               aux(i,j,2) = cl
            else
               aux(i,j,1) = rhor
               aux(i,j,2) = cr
            endif
         enddo
      enddo
      return
      end
