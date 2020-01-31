C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "CONSTANTS.H"
#include "REAL.H"

      subroutine LEGENDRECOEFFS(coeffs, l, m)
!     Retrieves polynomial coefficients of an associated Legendre function.
!
!     The (l, m) associated Legendre function at x is
!     P_l^m(x) = ((-1)^m / (2^l * l!)) * (1 - x^2)^(m/2)
!                * d^{l+m}/dx^{l+m} ((x^2 - 1)^l) .
!     This function returns the coefficients of the polynomial
!     ((-1)^m / (2^l * l!)) * d^{l+m}/dx^{l+m} ((x^2 - 1)^l)
!     so that in order to get P_l^m(x), we evaluate this polynomial
!     and multiply by (1 - x^2)^(m/2); this is done in LEGENDREEVAL.
!     In order to get P_l^m'(x), we call LEGENDREDERIV.
!
!     => inputs,  <= outputs
!     coeffs <= Legendre coefficients, indexed by 0:l-m
!     l => positive integer index
!     m => index in range 0:l
!
!     petermc, 21 Feb 2008
!     --- arguments ---
      implicit none
      integer l, m
      REAL_T coeffs(0:l-m)
!     --- locals ---
      integer cpoly(0:2*l)
      integer degc, iter, i
      REAL_T fac

      if (.not. (l .ge. 1)) then
         print *, 'LEGENDRECOEFFS given l = ', l
         print *, 'LEGENDRECOEFFS requires l >= 1'
         call MAYDAY_ERROR()
      elseif (.not. ((m .ge. 0) .and. (m .le. l))) then
         print *, 'LEGENDRECOEFFS given l = ', l, ' m = ', m
         print *, 'LEGENDRECOEFFS requires 0 <= m <= l'
         call MAYDAY_ERROR()
      endif

!
!     First, fill cpoly(0:2*l) with 
!     the coefficients of the polynomial (-1 + x^2)^l.
!     (Recall l >= 1.)
!

!     Initialize with polynomial = -1 + x^2, degree 2.
      cpoly(0) = -1
      cpoly(1) = 0
      cpoly(2) = 1
      degc = 2
!     Multiply cpoly by (-1 + x^2), l-1 times
      do iter = 2, l
         cpoly(degc + 2) = cpoly(degc)
         cpoly(degc + 1) = cpoly(degc - 1)
         do i = degc, 2, -1
            cpoly(i) = -cpoly(i) + cpoly(i - 2)
         enddo
         cpoly(1) = -cpoly(1)
         cpoly(0) = -cpoly(0)
         degc = degc + 2
      enddo
!     Now cpoly(0:degc) contains coefficients of (-1 + x^2)^l,
!     and degc = 2*l.

!
!     Now replace cpoly(0:degc) by (l + m)th derivative of (-1 + x^2)^l.
!
      do iter = 1, l + m
!     Replace cpoly(0:degc) by its derivative.
         do i = 1, degc
            cpoly(i - 1) = i * cpoly(i)
         enddo
         degc = degc - 1
      enddo

!
!     Set fac = (-1)^m / (2^l * l!)
!

      fac = one
      if (mod(m, 2) .eq. 1) then
         fac = -one
      endif
      do i = 1, l
         fac = fac / (i * two)
      enddo

      do i = 0, degc
         coeffs(i) = cpoly(i) * fac
      enddo

!     Now to evaluate, just need to sum up coeffs * powers of x,
!     and multiply the result by (1 - x^2)^(m/2).
!     We do this in function LEGENDREEVAL.

      return
      end


      REAL_T function LEGENDREEVAL(x, coeffs, l, m)
!     Evaluates associated Legendre function at x:
!     P_l^m(x) = ((-1)^m / (2^l * l!)) * (1 - x^2)^(m/2)
!                * d^{l+m}/dx^{l+m} ((x^2 - 1)^l)
!              = (1 - x^2)^(m/2) * sum_{i = 0:l-m} (coeffs[i] * x^i)
!
!     => inputs,  <= outputs
!     coeffs => Legendre coefficients, indexed by 0:l-m, from LEGENDRECOEFFS
!     l => positive integer index
!     m => index in range 0:l
!
!     petermc, 21 Feb 2008
!     --- arguments ---
      implicit none
      integer l, m
      REAL_T x, coeffs(0:l-m)
!     --- locals ---
      integer i
      REAL_T ans

!     polynomial evaluation
      ans = coeffs(l - m)
      do i = l - m - 1, 0, -1
         ans = ans * x + coeffs(i)
      enddo
!     multiplication by (1 - x^2)^(m/2)
      if (mod(m, 2) .eq. 0) then
         ans = ans * (one - x*x)**(m/2)
      else
         ans = ans * (sqrt(one - x*x))**m
      endif

      legendreeval = ans

      return
      end


      REAL_T function LEGENDREDERIV(x, coeffs, l, m)
!     Evaluates derivative of associated Legendre function at x, which is:
!     P_l^m(x) = ((-1)^m / (2^l * l!)) * (1 - x^2)^(m/2)
!                * d^{l+m}/dx^{l+m} ((x^2 - 1)^l)
!              = (1 - x^2)^(m/2) * sum_{i = 0:l-m} (coeffs[i] * x^i)
!     and the derivative is:
!     P_l^m'(x) = ((-1)^m / (2^l * l!)) * (1 - x^2)^(m/2)
!                 * d^{l+m}/dx^{l+m} ((x^2 - 1)^l)
!
!     => inputs,  <= outputs
!     coeffs => Legendre coefficients, indexed by 0:l-m, from LEGENDRECOEFFS
!     l => positive integer index
!     m => index in range 0:l
!
!     petermc, 21 Feb 2008
!     --- arguments ---
      implicit none
      integer l, m
      REAL_T x, coeffs(0:l-m)
!     --- locals ---
      integer i
      REAL_T ans

!     polynomial evaluation
      ans = -l * coeffs(l - m)
      ans = ans * x + (1 - l) * coeffs(l - m - 1)
      do i = l - m - 1, 2, -1
         ans = ans * x + (1-m-i) * coeffs(i-1) + (i+1) * coeffs(i+1)
      enddo
      ans = ans * x + (2*coeffs(2) - m * coeffs(0))
      ans = ans * x + coeffs(1)
!     multiplication by (1 - x^2)^(m/2 - 1)
      if (mod(m, 2) .eq. 0) then
         ans = ans * (one - x*x)**(m/2 - 1)
      else
         ans = ans * (sqrt(one - x*x))**(m - 2)
      endif

      legendrederiv = ans

      return
      end
