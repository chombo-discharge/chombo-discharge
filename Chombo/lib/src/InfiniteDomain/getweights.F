!    _______              __
!   / ___/ /  ___  __ _  / /  ___
!  / /__/ _ \/ _ \/  V \/ _ \/ _ \
!  \___/_//_/\___/_/_/_/_.__/\___/
!
!  CHOMBO Copyright (c) 2000-2004, The Regents of the University of
!  California, through Lawrence Berkeley National Laboratory (subject to
!  receipt of any required approvals from U.S. Dept. of Energy).  All
!  rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!
!  (1) Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!  (2) Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!  (3) Neither the name of Lawrence Berkeley National Laboratory, U.S.
!  Dept. of Energy nor the names of its contributors may be used to endorse
!  or promote products derived from this software without specific prior
!  written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!  You are under no obligation whatsoever to provide any bug fixes,
!  patches, or upgrades to the features, functionality or performance of
!  the source code ("Enhancements") to anyone; however, if you choose to
!  make your Enhancements available either publicly, or directly to
!  Lawrence Berkeley National Laboratory, without imposing a separate
!  written license agreement for such Enhancements, then you hereby grant
!  the following license: a non-exclusive, royalty-free perpetual license
!  to install, use, modify, prepare derivative works, incorporate into
!  other computer software, distribute, and sublicense such Enhancements or
!  derivative works thereof, in binary and source code form.
!
!  TRADEMARKS. Product and company names mentioned herein may be the
!  trademarks of their respective owners.  Any rights not expressly granted
!  herein are reserved.

#include "CONSTANTS.H"

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
      implicit none
!     --- arguments ---
      integer lo, hi, intrul, ier
      real*8 wtint(lo:hi), h
!     --- locals ---
      real*8 fac
      integer k, ptchsz

      ptchsz = hi - lo

!     wtint:  integration weights

!     if (ptchsz .eq. 0) then perhaps do same as if (intrul .eq. 0)

      if (intrul .eq. 0) then
         do k = lo, hi
            wtint(k) = 1d0
         enddo

      elseif (intrul .eq. 1) then
!     trapezoidal rule:  1, 1: * 1/2

         wtint(lo) = h / 2d0
         do k = lo + 1, hi - 1
            wtint(k) = h
         enddo
         wtint(hi) = h / 2d0

      elseif (intrul .eq. 2) then
!     Simpson's rule:  1, 4, 1: * 1/3
         if (mod(ptchsz, 2) .ne. 0) then
            print *, 'getintegral1dweights:  Simpson rule (2) ',
     &           'has length ', ptchsz, ' indivisible by 2'
            ier = -2
            return
         endif

         fac = h / 3d0
         wtint(lo) = fac
         do k = lo + 1, hi - 1, 2
            wtint(k) = 4d0 * fac
         enddo
         do k = lo + 2, hi - 2, 2
            wtint(k) = 2d0 * fac
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

         fac = h * (2d0/45d0)
         wtint(lo) = 7d0 * fac
         do k = lo + 1, hi - 1, 2
            wtint(k) = 32d0 * fac
         enddo
         do k = lo + 2, hi - 2, 4
            wtint(k) = 12d0 * fac
         enddo
         do k = lo + 4, hi - 4, 4
            wtint(k) = 14d0 * fac
         enddo
         wtint(hi) = 7d0 * fac

      elseif (intrul .eq. 8) then
!     Newton-Cotes for 9 points:
!     989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989: * 4/14175
         if (mod(ptchsz, 8) .ne. 0) then
            print *, 'getintegral1dweights:  Newton-Cotes rule (8) ',
     &           'has length ', ptchsz, ' indivisible by 8'
            ier = -8
            return
         endif

         fac = h * (4d0/14175d0)
         wtint(lo) = 989d0 * fac
         do k = lo + 1, hi - 7, 8
            wtint(k) = 5888d0 * fac
         enddo
         do k = lo + 7, hi - 1, 8
            wtint(k) = 5888d0 * fac
         enddo
         do k = lo + 2, hi - 6, 8
            wtint(k) = -928d0 * fac
         enddo
         do k = lo + 6, hi - 2, 8
            wtint(k) = -928d0 * fac
         enddo
         do k = lo + 3, hi - 5, 8
            wtint(k) = 10496d0 * fac
         enddo
         do k = lo + 5, hi - 3, 8
            wtint(k) = 10496d0 * fac
         enddo
         do k = lo + 4, hi - 4, 8
            wtint(k) = -4540d0 * fac
         enddo
         do k = lo + 8, hi - 8, 8
            wtint(k) = 2d0 * 989d0 * fac
         enddo
         wtint(hi) = 989d0 * fac

      else
         print *, 'getintegral1dweights:  illegal integration rule ',
     &        intrul
         ier = -99
         return
      endif

      ier = 0

      return
      end
