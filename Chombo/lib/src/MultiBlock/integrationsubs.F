C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "CONSTANTS.H"
#include "REAL.H"

      subroutine getintegral1dweights(wtint, lo, hi, h, intrul, ier)
!     This subroutine returns the weights for a one-dimensional integral.
!
!     => inputs, <= outputs
!     wtint <= array, weights
!     lo => index of low point of wtint
!     hi => index of high point of wtint
!     h => mesh spacing
!     intrul => which integration rule to use.
!              0 for constant;
!              1 for trapezoidal;
!              2 for Simpson, requiring hi-lo divisible by 2;
!              4 for Boole, requiring hi-lo divisible by 4;
!              8 for Newton-Cotes-9, requiring hi-lo divisible by 8.
!     ier <= error code.
!           0 if successful;
!           -2 if Simpson's rule and hi-lo not divisible by 2;
!           -4 if Boole's rule and hi-lo not divisible by 4;
!           -8 if Newton-Cotes-9 rule and hi-lo not divisible by 8;
!           -99 if illegal intrul
      implicit none
!     --- arguments ---
      integer lo, hi, intrul, ier
      REAL_T wtint(lo:hi), h
!     --- locals ---
      REAL_T fac
      integer k, ptchsz

      ptchsz = hi - lo

!     wtint:  integration weights

!     if (ptchsz .eq. 0) then perhaps do same as if (intrul .eq. 0)

      if (intrul .eq. 0) then
         do k = lo, hi
            wtint(k) = one
         enddo

      elseif (intrul .eq. 1) then
!     trapezoidal rule:  1, 1: * 1/2

         wtint(lo) = half * h
         do k = lo + 1, hi - 1
            wtint(k) = h
         enddo
         wtint(hi) = half * h

      elseif (intrul .eq. 2) then
!     Simpson's rule:  1, 4, 1: * 1/3
         if (mod(ptchsz, 2) .ne. 0) then
            print *, 'getintegral1dweights:  Simpson rule (2) ',
     &           'has length ', ptchsz, ' indivisible by 2'
            ier = -2
            return
         endif

         fac = third * h
         wtint(lo) = fac
         do k = lo + 1, hi - 1, 2
            wtint(k) = four * fac
         enddo
         do k = lo + 2, hi - 2, 2
            wtint(k) = two * fac
         enddo
         wtint(hi) = fac

      elseif (intrul .eq. 4) then
!     Boole's rule:  7, 32, 12, 32, 7: * 2/45
         if (mod(ptchsz, 4) .ne. 0) then
            print *, 'getintegral1dweights:  Boole rule (4) ',
     &           'has length ', ptchsz, ' indivisible by 4'
            ier = -4
            return
         endif

         fac = h * (two / (45*one))
         wtint(lo) = 7 * fac
         do k = lo + 1, hi - 1, 2
            wtint(k) = 32 * fac
         enddo
         do k = lo + 2, hi - 2, 4
            wtint(k) = 12 * fac
         enddo
         do k = lo + 4, hi - 4, 4
            wtint(k) = 14 * fac
         enddo
         wtint(hi) = 7 * fac

      elseif (intrul .eq. 8) then
!     Newton-Cotes for 9 points:
!     989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989: * 4/14175
         if (mod(ptchsz, 8) .ne. 0) then
            print *, 'getintegral1dweights:  Newton-Cotes rule (8) ',
     &           'has length ', ptchsz, ' indivisible by 8'
            ier = -8
            return
         endif

         fac = h * (four/(14175*one))
         wtint(lo) = 989 * fac
         do k = lo + 1, hi - 7, 8
            wtint(k) = 5888 * fac
         enddo
         do k = lo + 7, hi - 1, 8
            wtint(k) = 5888 * fac
         enddo
         do k = lo + 2, hi - 6, 8
            wtint(k) = -928 * fac
         enddo
         do k = lo + 6, hi - 2, 8
            wtint(k) = -928 * fac
         enddo
         do k = lo + 3, hi - 5, 8
            wtint(k) = 10496 * fac
         enddo
         do k = lo + 5, hi - 3, 8
            wtint(k) = 10496 * fac
         enddo
         do k = lo + 4, hi - 4, 8
            wtint(k) = -4540 * fac
         enddo
         do k = lo + 8, hi - 8, 8
            wtint(k) = two * 989 * fac
         enddo
         wtint(hi) = 989 * fac

      else
         print *, 'getintegral1dweights:  illegal integration rule ',
     &        intrul
         ier = -99
         return
      endif

      ier = 0

      return
      end


      subroutine getgaussquad1d(points, wts, ilo, ihi, alo, ahi)
!     This subroutine returns the points and weights for Gaussian
!     quadrature in one dimension.
!
!     => inputs, <= outputs
!     points <= array, evaluation points
!     wtint <= array, weights at evaluation points
!     ilo => low index of points and wtint
!     ihi => high index of points and wtint
!     alo => lower limit of interval of integration
!     ahi => upper limit of interval of integration
      implicit none
!     --- arguments ---
      integer ilo, ihi
      REAL_T points(ilo:ihi), wts(ilo:ihi), alo, ahi
!     --- locals ---
      REAL_T r1, r2, r3, r4, w1, w2, w3, w4, alenhf, amid
      integer npts, ipt

      npts = ihi - ilo + 1

!     First set points in [-1:1].
      if (npts .eq. 1) then
         points(ilo) = zero
         wts(ilo) = two
      elseif (npts .eq. 2) then
         r1 = one / sqrt(three)
         points(ilo) = -r1
         points(ilo+1) = r1
         wts(ilo) = one
         wts(ilo+1) = one
      elseif (npts .eq. 3) then
         r1 = sqrt(three / five)
         points(ilo) = -r1
         points(ilo+1) = zero
         points(ilo+2) = r1
         w1 = eight / nine
         w2 = five /nine
         wts(ilo) = w2
         wts(ilo+1) = w1
         wts(ilo+2) = w2
      elseif (npts .eq. 4) then
         r1 = three / seven
         r2 = two * sqrt(six / five) / seven
         r3 = sqrt(r1 - r2)
         r4 = sqrt(r1 + r2)
         points(ilo) = -r4
         points(ilo+1) = -r3
         points(ilo+2) = r3
         points(ilo+3) = r4
         w1 = half
         w2 = sqrt(30*one) / (36*one)
         w3 = w1 + w2
         w4 = w1 - w2
         wts(ilo) = w4
         wts(ilo+1) = w3
         wts(ilo+2) = w3
         wts(ilo+3) = w4
      elseif (npts .eq. 5) then
         r1 = five
         r2 = two * sqrt(ten / seven)
         r3 = sqrt(r1 - r2) * third
         r4 = sqrt(r1 + r2) * third
         points(ilo) = -r4
         points(ilo+1) = -r3
         points(ilo+2) = zero
         points(ilo+3) = r3
         points(ilo+4) = r4
         w1 = (322*one) / (900*one)
         w2 = (13*one) * sqrt(70*one) / (900*one)
         w3 = w1 + w2
         w4 = w1 - w2
         wts(ilo) = w4
         wts(ilo+1) = w3
         wts(ilo+2) = (128*one) / (225*one)
         wts(ilo+3) = w3
         wts(ilo+4) = w4
      else
         print *, 'GETGAUSSQUAD1D points only 1 through 5, not ',
     &        npts
         call MAYDAY_ERROR()
      endif

!     Now change interval from [-1:1] to [alo:ahi].
      alenhf = (ahi - alo) * half
      amid = (alo + ahi) * half
      do ipt = ilo, ihi
         points(ipt) = amid + points(ipt) * alenhf
         wts(ipt) = alenhf * wts(ipt)
      enddo

      return
      end
