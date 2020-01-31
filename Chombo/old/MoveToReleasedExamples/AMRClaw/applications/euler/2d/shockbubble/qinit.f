c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit none

       integer maxmx,maxmy,meqn, mbc,mx,my, maux
       double precision xlower, ylower,dx,dy

       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       double precision qin(5),qout(5), rinf, vinf, einf
       double precision xlow, ylow, win
       integer i,j, m

       common /comic/ qin, qout
       common /cominf/ rinf,vinf,einf

       do i = 1,mx
          xlow = xlower + (i-1)*dx
          do j = 1,my
c            # set (xlow,ylow) to lower left corner of grid cell:
             ylow = ylower + (j-1)*dy
             call cellave2(xlow,ylow,dx,dy,win)
c            # win is now the fraction of the cell that lies inside the circle
             do  m = 1,meqn
                q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
             enddo
          enddo
          if (xlow .lt. 0.2d0) then
c            # behind shock:
             do j = 1,my
                q(i,j,1) = rinf
                q(i,j,2) = rinf*vinf
                q(i,j,3) = 0.d0
                q(i,j,4) = einf
                q(i,j,5) = 0.d0
             enddo
          endif
       enddo

       return
       end
