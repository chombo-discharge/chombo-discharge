C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "REAL.H"
#include "CONSTANTS.H"

      subroutine equiangulartolonlat(
     &     lon, lat, xi, eta, panel, eps)
!     Computes longitude and latitude from mapped coordinates
!
!     => inputs,  <= outputs
!     lon <= longitude, in radians
!     lat <= latitude, in radians
!     xi => mapped first coordinate
!     eta => mapped second coordinate
!     panel => panel number, 0 thru 5
!     eps => machine precision, DBL_EPSILON
      implicit none
!     --- arguments ---
      real_t lon, lat, xi, eta, eps
      integer panel
!     --- locals ---
      real_t dgx, dgy

      dgx = tan(xi - fourth*Pi)
      dgy = tan(eta - fourth*Pi)

      if ((panel .ge. 0) .and. (panel .le. 3)) then
!     equatorial panels
         lon = atan(dgx) + panel * half*Pi
         lat = atan(dgy / sqrt(one + dgx*dgx))
      elseif (panel .eq. 4) then
!     north polar panel
         if ((abs(dgx) .gt. eps) .and. (abs(dgy) .gt. eps)) then
            lon = atan2(dgx, -dgy)
         elseif ((abs(dgx) .le. eps) .and. (dgy .ge. eps)) then
            lon = Pi
         elseif ((abs(dgy) .le. eps) .and. (dgx .ge. eps)) then
            lon = half*Pi
         elseif ((abs(dgy) .le. eps) .and. (dgx .le. -eps)) then
            lon = three*half*Pi
         else
            lon = zero
         endif
         lat = half*Pi - atan(sqrt(dgx*dgx + dgy*dgy))
      elseif (panel .eq. 5) then
!     south polar panel
         if ((abs(dgx) .gt. eps) .and. (abs(dgy) .gt. eps)) then
            lon = atan2(dgx, dgy)
         elseif ((abs(dgx) .le. eps) .and. (dgy .le. -eps)) then
            lon = Pi
         elseif ((abs(dgy) .le. eps) .and. (dgx .ge. eps)) then
            lon = half*Pi
         elseif ((abs(dgy) .le. eps) .and. (dgx .le. -eps)) then
            lon = three*half*Pi
         else
            lon = zero
         endif
         lat = -half*Pi + atan(sqrt(dgx*dgx + dgy*dgy))
      else
         print *, 'Error: invalid panel ', panel
         call MAYDAY_ERROR()
      endif

      return
      end


      subroutine equiangulartocartesian(
     &     x, y, z, xi, eta, panel)
!     Computes Cartesian coordinates from mapped coordinates
!
!     => inputs,  <= outputs
!     x <= cartesian x
!     y <= cartesian y
!     z <= cartesian z
!     xi => mapped first coordinate
!     eta => mapped second coordinate
!     panel => panel number, 0 thru 5
      implicit none
!     --- arguments ---
      real_t x, y, z, xi, eta
      integer panel
!     --- locals ---
      real_t dgx, dgy, dsx, dsy, dsz

      dgx = tan(xi - fourth*Pi)
      dgy = tan(eta - fourth*Pi)

      dsz = one / sqrt(one + dgx*dgx + dgy*dgy)
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


      subroutine vectorcartesiantoequiangular(
     &     tvxi, tveta, xi, eta, panel, vx, vy, vz)
!     Transforms vector from real-coordinate basis to mapped-coordinate basis
!
!     => inputs,  <= outputs
!     tvxi, tveta <= components of vector in mapped-coordinate basis
!     xi, eta => position in mapped coordinates
!     panel => panel number, 0 thru 5
!     vx, vy, vz => components of vector in real-coordinate basis
      implicit none
!     --- arguments ---
      real_t tvxi, tveta, xi, eta, vx, vy, vz
      integer panel
!     --- locals ---
      real_t dgx, dgy, delta, denomx, denomy, onex, oney, fracx, fracy

      dgx = tan(xi - fourth*Pi)
      dgy = tan(eta - fourth*Pi)

      delta = sqrt(one + dgx*dgx + dgy*dgy)

      onex = delta / (one + dgx*dgx)
      oney = delta / (one + dgy*dgy)

      fracx = dgx * onex
      fracy = dgy * oney

      if (panel .eq. 0) then
         tvxi = -fracx * vx + onex * vy
         tveta = -fracy * vx + oney * vz
      elseif (panel .eq. 1) then
         tvxi = -onex * vx - fracx * vy
         tveta = -fracy * vy + oney * vz
      elseif (panel .eq. 2) then
         tvxi = fracx * vx - onex * vy
         tveta = fracy * vx + oney * vz
      elseif (panel .eq. 3) then
         tvxi = onex * vx + fracx * vy
         tveta = fracy * vy + oney * vz
      elseif (panel .eq. 4) then
         tvxi = onex * vy - fracx * vz
         tveta = -oney * vx - fracy * vz
      elseif (panel .eq. 5) then
         tvxi = onex * vy + fracx * vz
         tveta = oney * vx + fracy * vz
      else
         print *, 'Error: invalid panel ', panel
         call MAYDAY_ERROR()
      endif

      return
      end
