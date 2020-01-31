      subroutine mapc2p(xc,yc,xp,yp)
      implicit none

      double precision xc,yc,xp,yp
      double precision pi, theta, r

      common /compi/ pi

      theta = 2.d0*pi*yc
      r = 16.d0*xc + 1.d0

      xp = r * cos(theta)
      yp = r * sin(theta)

      return
      end
