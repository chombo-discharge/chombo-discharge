c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &                   xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit none
c

       integer i,j,k,mbc,maxmx, maxmy, mx, my, maux, meqn
       double precision xlower, ylower, dx,dy
       double precision gamma, gamma1,ur, ul, rhor, rhol,
     &       pr, pl, rhour, rhoul, el, er, win

       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

       double precision xe,xc,ye,yc,xp,yp

       common /param/ gamma, gamma1
       common /indata/ rhol,rhor,pl,pr,ul,ur

       do j = 1 - mbc, my + mbc
          ye = ylower + (j-1.d0)*dy
          yc = ylower + (j-0.5d0)*dy
          do i = 1 - mbc, mx + mbc
             xe = xlower + (i-1.d0)*dx
             xc = xlower + (i-0.5d0)*dx

c            # cellave2 calls fdisc (below).
c            # fdisc < 0 --> win = 1.0
c            # fdisc > 0 --> win = 0.0
             call cellave2(xe,ye,dx,dy,win)

             el = pl/gamma1 + 0.5d0*rhol*ul*ul
             er = pr/gamma1 + 0.5d0*rhor*ur*ur

             rhoul = rhol*ul
             rhour = rhor*ur

             q(i,j,1) = win*rhol + (1.d0-win)*rhor
             q(i,j,2) = win*rhoul + (1.d0-win)*rhour
             q(i,j,3) = 0.d0
             q(i,j,4) = win*el + (1.d0-win)*er
          enddo
       enddo

       return
       end


      double precision function fdisc(x,y)
      implicit none

      double precision x,y,xp,yp

      call mapc2p(x,y,xp,yp)

c     # x < -2 is negative; area fraction = 1
      fdisc = xp + 2.0d0

      return
      end
