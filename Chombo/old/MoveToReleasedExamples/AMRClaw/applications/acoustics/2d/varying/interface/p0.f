c     =====================================================
      double precision function p0(r)
c     =====================================================
c
      implicit none

      double precision width, s, pi, r

      common /compi/ pi

      width = 0.1d0
      s = 0.d0
      if (abs(r-s) .lt. width) then
              p0 = 1.d0 + 0.5d0 * (dcos(pi*(r-s)/width) - 1.d0)
         else
              p0 = 0.d0
         endif
      return
      end
