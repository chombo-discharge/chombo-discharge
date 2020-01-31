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
      subroutine deriv1overr(cpind, cfac, cpow, nterms, deg)
!     Get expressions for partial derivatives of 1/(4*pi*r) wrt x and y.
!     For (d^(p1+p2) / dx^p1 dy^p2) (1/(4*pi*r)),
!     set pind = (p1*(2*deg+3-p1))/2 + p2 + 1
!     and then the expression for this (p1, p2) partial derivative is
!     sum_{d=0:nterms, cpind(d)==pind}
!          cfac(d) * (x^cpow(1, d)) * (y^cpow(2, d)) / (r^(1+2*cpow(3, d))).
!     where r = sqrt(x^2 + y^2).
!
!     nterms:  output, the number of terms
!     deg:  input, the max degree of partial derivatives to use; we take
!           d^(p1+p2) / dx^p1 dy^p2, for all p1 >= 0, p2 >= 0 provided
!           p1 + p2 <= deg
!     cpind:  output, index of the derivative for which the term applies
!     cfac:  output, the constant factor of the term
!     cpow:  output, the powers of x, y, r in the term
      implicit none
!     --- arguments ---
      integer nterms, deg
      integer cpind(((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
      integer cpow(3, ((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
      real*8 cfac(((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
!     --- locals ---
      integer p1, p2, pind, prvxf, prvxl, prvyf, prvyl, term
      real*8 pi, datan, fac

      pi = 4 * datan(1d0)
      fac = 1d0 / (4*pi)

!     Terms for a particular partial derivative are stored consecutively
!     in the cfac and cpow arrays:  thus, the cpind array has a bunch of
!     1's followed by a bunch of 2's followed by a bunch of 3's, etc.

!     Example:  deg == 4, the index pind = (p1*(2*deg+3-p1))/2 + p2 + 1
!     is shown in the table.
!
!      |  4 | 5
!      |  3 | 4  9
!     p2  2 | 3  8  12
!      |  1 | 2  7  11  14
!      |  0 | 1  6  10  13  15
!            -----------------
!             0  1   2   3   4
!            ------ p1 -------

!     prvxf, prvxl are indices of first and last terms for
!     previous partial wrt x only; that is, for
!     (d^(p1-1) / dx^(p1-1)) (1/(4*pi*r)).
!     Then we differentiate this expression wrt x to get
!     (d^p1 / dx^p1) (1/(4*pi*r)).

!     prvyf, prvyl are indices of first and last terms for
!     expression for (d^(p1+p2-1) / dx^p1 dy^(p2-1)) (1/(4*pi*r)).
!     Then we differentiate this expression wrt y to get
!     (d^(p1+p2) / dx^p1 dy^p2) (1/(4*pi*r)).

!     number of terms used in all expressions.
!     increment this when you get a new term.
      nterms = 0
!     pind = (p1*(2*deg+3-p1))/2 + p2 + 1
!     ranges from 1 to ((deg+1)*(deg+2))/2
      pind = 0
      do p1 = 0, deg
         do p2 = 0, deg
            if (p1 + p2 .le. deg) then
               pind = pind + 1
!     now find expression for (d^(p1+p2) / dx^p1 dy^p2) (1/(4*pi*r))
               if (p1 .eq. 0 .and. p2 .eq. 0) then
!     (p1, p2) = (0, 0) corresponds to 1/(4*pi*r)
                  cpind(1) = pind
                  cfac(1) = fac
                  cpow(1, 1) = 0
                  cpow(2, 1) = 0
!     At this stage, cpow(3, :) indicates the power of r.
!     At the end of this subroutine, we modify it.  See below.
                  cpow(3, 1) = -1
!     increment nterms.  update prvxf, prvxl, prvyf, prvyl.
                  prvxf = 1
                  prvxl = 1
                  prvyf = 1
                  prvyl = 1
                  nterms = nterms + 1
               else
!     get new expression by differentiating previous one.
                  if (p2 .eq. 0) then
!     next partial wrt x:
!     differentiate expression from prvxf to prvxl by x,
!     and update nterms, prvxf and prvxl while you're at it.
                     call diffexpr(cpind, cfac, cpow, deg, nterms,
     &                    prvxf, prvxl, 1, pind)
                     nterms = prvxl
!     set so that this is next expression differentiated by y
                     prvyf = prvxf
                     prvyl = prvxl
                  else
!     next partial wrt y:
!     differentiate expression from prvyf to prvyl by y,
!     and update nterms, prvyf and prvyl while you're at it.
                     call diffexpr(cpind, cfac, cpow, deg, nterms,
     &                    prvyf, prvyl, 2, pind)
                     nterms = prvyl
                  endif
               endif
            endif
         enddo
      enddo
!     Since cpow(3, :) must always be negative and odd, transform it.
!     Old:  r^cpow(3, :)
!     New:  1/r^(1 + 2*cpow(3, :))
      do term = 1, nterms
         cpow(3, term) = (-1 - cpow(3, term)) / 2
      enddo

      return
      end

c=======================================================================

      subroutine diffexpr(cpind, cfac, cpow, deg, nterms,
     &     prvf, prvl, ivar, pind)
!     This subroutine is called by deriv1overr.
!     It appends to the arrays cpind, cfac, cpow
!     with terms for next partial derivatives.
      implicit none
!     --- arguments ---
      integer deg, nterms, prvf, prvl, ivar, pind
      integer cpind(((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
      real*8 cfac(((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
      integer cpow(3, ((deg+1)*(deg+2)*(deg+3)*(deg+4))/24)
!     --- locals ---
      integer powcur(3), pownew(3), term, termo, ntold
      real*8 faccur, facnew
      logical found

      ntold = nterms
      do term = prvf, prvl
         faccur = cfac(term)
         powcur(1) = cpow(1, term)
         powcur(2) = cpow(2, term)
         powcur(3) = cpow(3, term)
!     differentiate wrt x:
!     d/dx (fac * x^pow(1) * y^pow(2) * r^pow(3)) =
!     fac * pow(1) * x^(pow(1)-1) * y^pow(2) * r^pow(3) +
!     fac * pow(3) * x^(pow(1)+1) * y^pow(2) * r^(pow(3)-2).
!     differentiate wrt y:
!     d/dy (fac * x^pow(1) * y^pow(2) * r^pow(3)) =
!     fac * pow(2) * x^pow(1) * y^(pow(2)-1) * r^pow(3) +
!     fac * pow(3) * x^pow(1) * y^(pow(2)+1) * r^(pow(3)-2).
         if (powcur(ivar) .ne. 0) then
!     first term
            facnew = powcur(ivar) * faccur
!     power of differentiating variable in the term
            pownew(ivar) = powcur(ivar) - 1
!     power of the other variable in the term, remains unchanged
            pownew(3-ivar) = powcur(3-ivar)
!     power of r in the term
            pownew(3) = powcur(3)
!     now incorporate this new term into the arrays.
!     if, in this call, we've already had a term with the same powers,
!     then change its coefficient.
            found = .false.
!     termo loops over indices of all terms added in this call
            do termo = ntold + 1, nterms
               if ((pownew(1) .eq. cpow(1, termo)) .and.
     &              (pownew(2) .eq. cpow(2, termo)) .and.
     &              (pownew(3) .eq. cpow(3, termo))) then
                  found = .true.
                  cfac(termo) = cfac(termo) + facnew
               endif
            enddo
            if (.not. found) then
!     new term, so need to append a new element to the arrays
               nterms = nterms + 1
               cpind(nterms) = pind
               cpow(1, nterms) = pownew(1)
               cpow(2, nterms) = pownew(2)
               cpow(3, nterms) = pownew(3)
               cfac(nterms) = facnew
            endif
         endif
         if (powcur(3) .ne. 0) then
!     second term
            facnew = powcur(3) * faccur
!     power of differentiating variable
            pownew(ivar) = powcur(ivar) + 1
!     power of the other variable, remains unchanged
            pownew(3-ivar) = powcur(3-ivar)
!     power of r
            pownew(3) = powcur(3) - 2
!     now incorporate this new term into the arrays.
!     if, in this call, we've already had a term with the same powers,
!     then change its coefficient.
            found = .false.
!     termo loops over indices of all terms added in this call
            do termo = ntold + 1, nterms
               if ((pownew(1) .eq. cpow(1, termo)) .and.
     &              (pownew(2) .eq. cpow(2, termo)) .and.
     &              (pownew(3) .eq. cpow(3, termo))) then
                  found = .true.
                  cfac(termo) = cfac(termo) + facnew
               endif
            enddo
            if (.not. found) then
!     new term, so need to append a new element to the arrays
               nterms = nterms + 1
               cpind(nterms) = pind
               cpow(1, nterms) = pownew(1)
               cpow(2, nterms) = pownew(2)
               cpow(3, nterms) = pownew(3)
               cfac(nterms) = facnew
            endif
         endif

      enddo

!     set prvf to index of first term of new expression
      prvf = ntold + 1
!     set prvl to index of last term of new expression
      prvl = nterms

      return
      end
