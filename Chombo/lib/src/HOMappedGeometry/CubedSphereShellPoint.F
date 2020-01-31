C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "REAL.H"
#include "CONSTANTS.H"

      real_t function shellface2dintegral(alpha, y)
!
!     This unlikely-looking function is the 2D integral of
!     (1 + X^2) * (1 + Y^2) / (1 + X^2 + Y^2)^(3/2),
!     where X = tan(alpha), Y = tan(beta).
!     Its derivative wrt alpha is Y / sqrt(1 + X^2 + Y^2),
!     and derivative wrt beta is X / sqrt(1 + X^2 + Y^2).
!     Found with Mathematica.
!
      implicit none
      real_t alpha, y, x, s2y1x, seca

      x = tan(alpha)
      seca = one / cos(alpha)
      s2y1x = sqrt(two * (one + y*y / (one + x*x)))
      shellface2dintegral =
     &     atan( (sqrt(two) * y * sin(alpha)) / s2y1x ) *
     &     s2y1x * seca / sqrt(two * (one + x*x + y*y))

      return
      end


      real_t function shellfaceintegral(alo, ahi, ylo, yhi)
      implicit none
      real_t alo, ahi, ylo, yhi
      real_t shellface2dintegral

      shellfaceintegral =
     &     (shellface2dintegral(ahi, yhi)
     &     - shellface2dintegral(ahi, ylo)) -
     &     (shellface2dintegral(alo, yhi)
     &     - shellface2dintegral(alo, ylo))

      return
      end

      subroutine equiangulartocartesian3d(
     &     x, y, z, xi, eta, r, panel, radius, height)
!     Computes Cartesian coordinates from mapped coordinates
!
!     => inputs,  <= outputs
!     x <= cartesian x
!     y <= cartesian y
!     z <= cartesian z
!     xi => mapped first coordinate
!     eta => mapped second coordinate
!     r => mapped third coordinate (note - no dilation with r)
!     panel => panel number, 0 thru 5
      implicit none
!     --- arguments ---
      real_t x, y, z, xi, eta, r, radius, height
      integer panel
!     --- locals ---
      real_t dgx, dgy, dsx, dsy, dsz, r0

      dgx = tan(xi - fourth*Pi)
      dgy = tan(eta - fourth*Pi)

!     Radius starts at 1 in cubed sphere shell model
!      r0 = 1
      r0=radius
      dsz = (r0 + height * r) / sqrt(one + dgx*dgx + dgy*dgy)
      dsx = dsz * dgx
      dsy = dsz * dgy

      if (panel .eq. 0) then
         x = dsz
         y = dsx
         z = dsy
      elseif (panel .eq. 1) then
         x = -dsx
         y = dsz
         z = dsy
      elseif (panel .eq. 2) then
         x = -dsz
         y = -dsx
         z = dsy
      elseif (panel .eq. 3) then
         x = dsx
         y = -dsz
         z = dsy
      elseif (panel .eq. 4) then
         x = -dsy
         y = dsx
         z = dsz
      elseif (panel .eq. 5) then
         x = dsy
         y = dsx
         z = -dsz
      else
         print *, 'Error: invalid panel ', panel
         call MAYDAY_ERROR()
      endif

      return
      end

