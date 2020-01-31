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
      subroutine coeffs1sidedderiv(coeff, d, h)
!     Returns coefficients of one-sided derivative of degree d,
!     by differentiating the Lagrange interpolating polynomial of
!     degree d.  The error is then O(h^d).
!
!     => input, <= outputs
!     coeff <= coefficients
!     d => degree of approximation, at least 1
!     h => mesh spacing
!
!     Then the formula for the derivative of f(a) is
!     sum_{i=0:d} coeff(i) * f(a + i*h)
!     sum_{i=0:d} -coeff(i) * f(a - i*h)
      implicit none
!     --- arguments ---
      integer d
      real*8 coeff(0:d)
      real*8 h
!     --- locals ---
      integer i, j, k
      real*8 sum, prod

      do i = 0, d
         sum = 0d0
         do j = 0, d
            if (j .ne. i) then
               prod = 1d0 / dfloat(i-j)
               do k = 0, d
                  if (k .ne. i .and. k .ne. j) then
                     prod = prod * (-dfloat(k) / dfloat(i-k))
                  endif
               enddo
               sum = sum + prod
            endif
         enddo
         coeff(i) = sum / h
      enddo

      return
      end


      subroutine coeffsderiv(coeff, d, h, start, evalpt)
!     Returns coefficients of derivative of degree d,
!     by differentiating the Lagrange interpolating polynomial of
!     degree d.  The error is then O(h^d).
!
!     => input, <= outputs
!     coeff <= coefficients
!     d => degree of approximation, at least 1
!     h => mesh spacing
!     start => where interpolating points are start + i*h, i=0:d
!     evalpt => where derivative is evaluated
!
!     Then the formula is
!     f'(evalpt) = sum_{i=0:d} coeff(i) * f(start + i*h)
!     f'(evalpt) = sum_{i=0:d} -coeff(i) * f(start - i*h)
      implicit none
!     --- arguments ---
      integer d
      real*8 coeff(0:d)
      real*8 h
      real*8 start, evalpt
!     --- locals ---
      integer i, j, k
      real*8 sum, prod, thispt

      do i = 0, d
         sum = zero
         do j = 0, d
            if (j .ne. i) then
               prod = one / ((i - j) * h)
               do k = 0, d
                  thispt = start + k*h
                  if ((k .ne. i) .and. (k .ne. j)) then
                     prod = prod * (evalpt - thispt) /
     &                    ((i - k) * h)
                  endif
               enddo
               sum = sum + prod
            endif
         enddo
         coeff(i) = sum
      enddo

      return
      end
