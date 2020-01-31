      subroutine estimate_error2(mx,my,mbc,meqn,xlower, ylower,
     &      dx,dy, t, level, isBoundary, tol, q, error_measure)
      implicit none

      integer mx,my,meqn, level, mbc
      double precision dx,dy,xlower, ylower
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error_measure(1:mx,1:my)

      double precision gx,gy,tol, t, divx, divy, err
      integer i,j, isBoundary(4), im, ip, jm, jp,
     &        imin, imax, jmin,jmax, m
      double precision xi,yj,zk


      imin = isBoundary(1)
      imax = mx + 1-isBoundary(2)

      jmin = isBoundary(3)
      jmax = my + 1 - isBoundary(4)

      do i = 1,mx
         ip = min(i+1,imax)
         im = max(i-1,imin)
         divx = real(ip - im)
         xi = xlower + (i-0.5d0)*dx
         do j = 1,my
            jp = min(j+1,jmax)
            jm = max(j-1,jmin)
            divy = real(jp - jm)
            yj = ylower + (j - 0.5d0)*dy

c           # This is what AMRClaw does
            err = 0.d0
            do m = 1,meqn
               gx = abs(q(ip,j,m) - q(im,j,m))
               gy = abs(q(i,jp,m) - q(i,jm,m))
               err = dmax1(gx,gy,err)
            enddo
c           # Always refine in center.
            if ((abs(xi - 2.5) .le. 0.5d0
     &            .and. abs(yj - 2.5d0) .le. 0.5d0) .and. t < 0.1) then
               err = 1.0
            else
               err = abs(q(i,j,1))
            endif
            error_measure(i,j) = err

c           # And this is what Chombo does
c           gx = (q(ip,j,1) - q(im,j,1))/divx
c           gy = (q(i,jp,1) - q(i,jm,1))/divy
c           error_measure(i,j) = sqrt(gx*gx + gy*gy)
         enddo
      enddo
      end
