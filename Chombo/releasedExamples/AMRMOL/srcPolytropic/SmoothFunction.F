#include "REAL.H"
#include "CONSTANTS.H"
      real_t function smoothened(x)
!
!     This function returns:
!     cos(pi/2*x)^6, if |x| < 1;
!     0, otherwise.
!
!     This is a function that peaks at 1 when x = 0,
!     is zero for |x| >= 1,
!     and has a continuous 5th derivative.
      implicit none
      real_t x

      if (abs(x) .ge. one) then
         smoothened = zero
      else
         smoothened = (cos(Pi * half * x))**6
      endif

      return
      end
