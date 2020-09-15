      integer function coarsendown(n, d)
!
!     Returns floor(n/d).
!
      implicit none
      integer n, d, q

      q = n / d
      if (q * d .gt. n) q = q - 1
      coarsendown = q

      return
      end


      integer function modnn(n, d)
!
!     Returns n mod d, in range 0:d-1.
!
      implicit none
      integer n, d, q, m
      integer coarsendown

      q = coarsendown(n, d)
      m = n - q*d
      modnn = m

      return
      end
