      subroutine mapc2p(xc,yc,zc,xp,yp,zp)
      implicit none

      double precision xc,yc,zc,xp,yp,zp
      double precision Rot(3,3), xv(3)
      double precision x,y,z,d, skewness, radius, th, pi
      double precision rc,a,c,s, xc_old
      integer i,j

      common /compi/ pi

c     # ------------------------------------
c     # Linear mapping
c      xp = 2.d0*xc - 3.d0  !! in [-3,3]
c      yp = 3.d0*yc - 1.d0  !! in [-1,1]
c      zp = 5.d0*zc - 1.d0  !!m in [-1,1]

c     # -----------------------------------------------------------
c     # Rectangular grid with rotated inner sphere of radius xlen.
c
c     # The mesh data inside the sphere will be rotated about the vector
c     # (x,y,z), an angle is pi*skewness

c     # This is mainly for Chomboclaw...
c      xc_old = xc
c      xc = 2.d0*xc - 3.d0

       x = 0.5d0
       y = 0.5d0
       z = 0.5d0

c     # Normalize (x,y,z) vector
      d = sqrt(x*x + y*y + z*z)
      x = x/d
      y = y/d
      z = z/d

      skewness = 0.5d0
      radius = 0.8d0
      th = pi*skewness

      rc = sqrt(xc**2 + yc**2 + zc**2)
      if (rc .le. radius) then
         a = (1.d0 + cos(pi*rc/radius))/2.d0;
      else
         a = 0.d0
      endif

      c = cos(a*th)
      s = sin(a*th)

      Rot(1,1) = x**2*(1.d0-c) + c
      Rot(1,2) = x*y *(1.d0-c) - z*s
      Rot(1,3) = x*z *(1.d0-c) + y*s
      Rot(2,1) = x*y *(1.d0-c) + z*s
      Rot(2,2) = y**2*(1.d0-c) + c
      Rot(2,3) = y*z *(1.d0-c) - x*s
      Rot(3,1) = x*z *(1.d0-c) - y*s
      Rot(3,2) = y*z *(1.d0-c) + x*s
      Rot(3,3) = z**2*(1.d0-c) + c

      do i = 1,3
         xv(i) = Rot(i,1)*xc + Rot(i,2)*yc + Rot(i,3)*zc
      enddo
      xp = xv(1)
      yp = xv(2)
      zp = xv(3)

c      xc = xc_old

      end
