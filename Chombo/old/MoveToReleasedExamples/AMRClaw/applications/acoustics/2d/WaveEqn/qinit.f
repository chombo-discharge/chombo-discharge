c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Acoustics with smooth radially symmetric profile to test accuracy
c
       implicit none

       integer maxmx, maxmy, meqn, mbc, mx, my, maux
       double precision xlower,ylower, dx, dy
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       integer i,j
       double precision xcell, ycell, r, pressure, pi, xc, yc, width

       common /initq/ width
       common /compi/ pi

       xc = 2.5d0
       yc = 2.5d0

       do i = 1,mx
          xcell = xlower + (i - 0.5)*dx
          do j = 1,my
             ycell = ylower + (j - 0.5)*dy
             r = sqrt((xcell - xc)**2 + (ycell - yc)**2)

c             pressure = exp(-r**2/width**2)/(width**2)
c             phi = exp(-r**2/width**2)/(width**2)
c             u = phi*(-2*r/(width**2))*(0.5d0/2)*(2*(xcell-xc))
c             v = phi*(-2*r/(width**2))*(0.5d0/2)*(2*(ycell-yc))
c             pressure = 0.d0

             if (r .lt. width) then
                pressure = 1.d0
             else
                pressure = 0.d0
             endif
             q(i,j,1) = pressure
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
          enddo
       enddo

       end
