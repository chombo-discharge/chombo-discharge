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
      subroutine diffphi(phiex, phi, n, layer, h, vrbos, difnrm)
!     Write out differences between phiex and phi at interior points of
!     their domain, D2.
      implicit none
!     --- arguments ---
      integer n, layer
      real*8 phiex(-layer:n+layer, -layer:n+layer, -layer:n+layer)
      real*8 phi(-layer:n+layer, -layer:n+layer, -layer:n+layer)
      real*8 h
      integer vrbos
      real*8 difnrm(0:2)
!     --- locals ---
      integer i, j, k, lovld, hivld
      real*8 exhr, calchr, dif, difmax, exmax, clcmax
      real*8 exabs, clcabs, difabs
      real*8 sum1d, sum2d, sum1c, sum2c, sum1e, sum2e, dv
      real*8 norm1d, norm2d, norm1c, norm2c, norm1e, norm2e
      integer difi, difj, difk, exi, exj, exk, clci, clcj, clck
      real*8 dabs, dsqrt
      logical output

      output = .false.

!     low and high valid
      lovld = 1-layer
      hivld = n+layer-1
      if (output) then
         open (unit=2, file='outfile', status='new')
         write (2,*) lovld, hivld
      endif
!     difmax will be set to max(abs(phiex - phi)), and
!     max location at (difi, difj, difk).
      difmax = 0d0
      difi = 0
      difj = 0
      difk = 0
!     exmax will be set to max(abs(phiex)), and
!     max location at (exi, exj, exk).
      exmax = 0d0
      exi = 0
      exj = 0
      exk = 0
!     clcmax will be set to max(abs(phi)), and
!     max location at (clci, clcj, clck).
      clcmax = 0d0
      clci = 0
      clcj = 0
      clck = 0
!     sum1e, sum2e will be sums of phiex and phiex^2 over all points
      sum1e = 0d0
      sum2e = 0d0
!     sum1c, sum2c will be sums of phi and phi^2 over all points
      sum1c = 0d0
      sum2c = 0d0
!     sum1d, sum2d will be sums of abs(phiex-phi) and (phiex-phi)^2 over all points
      sum1d = 0d0
      sum2d = 0d0
      do i = lovld, hivld
         do j = lovld, hivld
            do k = lovld, hivld
               exhr = phiex(i, j, k)
               calchr = phi(i, j, k)
               dif = calchr - exhr
               difabs = dabs(dif)
               clcabs = dabs(calchr)
               exabs = dabs(exhr)
               sum1e = sum1e + exabs
               sum1c = sum1c + clcabs
               sum1d = sum1d + difabs
               sum2e = sum2e + (exabs*exabs)
               sum2c = sum2c + (clcabs*clcabs)
               sum2d = sum2d + (difabs*difabs)
               if (difabs .gt. difmax) then
                  difmax = difabs
                  difi = i
                  difj = j
                  difk = k
               endif
               if (exabs .gt. exmax) then
                  exmax = exabs
                  exi = i
                  exj = j
                  exk = k
               endif
               if (clcabs .gt. clcmax) then
                  clcmax = clcabs
                  clci = i
                  clcj = j
                  clck = k
               endif
               if (output) then
                  write (2,*) i, j, k, exhr, calchr, calchr-exhr
               endif
            enddo
         enddo
      enddo
      dv = h**3
!     find 1-norms
      norm1e = sum1e * dv
      norm1c = sum1c * dv
      norm1d = sum1d * dv
!     find 2-norms
      norm2e = dsqrt(sum2e * dv)
      norm2c = dsqrt(sum2c * dv)
      norm2d = dsqrt(sum2d * dv)
      if (output) then
         close (unit=2)
      endif
      if (vrbos .gt. 0) then
         print *, 'NORMS of exact solution, calculated solution, and ',
     &        'difference'
         print *, 'displaying:  ',
     &        'point of maximum, max norm, L1 norm, L2 norm'
         print *, 'exact norms:  (', exi, exj, exk, ') ',
     &        exmax, norm1e, norm2e
         print *, 'calc norms:  (', clci, clcj, clck, ') ',
     &        clcmax, norm1c, norm2c
         print *, 'diff norms:  (', difi, difj, difk, ') ',
     &        difmax, norm1d, norm2d
      endif
      difnrm(0) = difmax
      difnrm(1) = norm1d
      difnrm(2) = norm2d

      return
      end

c=======================================================================

      subroutine diffcrsephi(phiex, phi, nbigc, cl, h, d2lo, d2hi,
     &     vrbos, difnrm)
!     Write out differences between phiex and phi at points of their
!     domain, outside (d2lo:d2hi)^3
      implicit none
!     --- arguments ---
      integer nbigc, cl, d2lo, d2hi
      real*8 phiex(-cl:nbigc+cl, -cl:nbigc+cl, -cl:nbigc+cl)
      real*8 phi(-cl:nbigc+cl, -cl:nbigc+cl, -cl:nbigc+cl)
      real*8 h
      integer vrbos
      real*8 difnrm(0:2)
!     --- locals ---
      integer i, j, k, lovld, hivld
      real*8 exhr, calchr, dif, difmax, exmax, clcmax
      real*8 exabs, clcabs, difabs
      real*8 sum1d, sum2d, sum1c, sum2c, sum1e, sum2e, dv
      real*8 norm1d, norm2d, norm1c, norm2c, norm1e, norm2e
      integer difi, difj, difk, exi, exj, exk, clci, clcj, clck
      real*8 dabs, dsqrt
      logical output

      output = .false.

!     low and high valid
      lovld = -cl
      hivld = nbigc+cl
      if (output) then
         open (unit=2, file='outfilecoarse', status='new')
         write (2,*) lovld, hivld
      endif
!     difmax will be set to max(abs(phiex - phi)), and
!     max location at (difi, difj, difk).
      difmax = 0d0
      difi = 0
      difj = 0
      difk = 0
!     exmax will be set to max(abs(phiex)), and
!     max location at (exi, exj, exk).
      exmax = 0d0
      exi = 0
      exj = 0
      exk = 0
!     clcmax will be set to max(abs(phi)), and
!     max location at (clci, clcj, clck).
      clcmax = 0d0
      clci = 0
      clcj = 0
      clck = 0
!     sum1e, sum2e will be sums of phiex and phiex^2 over all points
      sum1e = 0d0
      sum2e = 0d0
!     sum1c, sum2c will be sums of phi and phi^2 over all points
      sum1c = 0d0
      sum2c = 0d0
!     sum1d, sum2d will be sums of abs(phiex-phi) and (phiex-phi)^2 over all points
      sum1d = 0d0
      sum2d = 0d0
      do i = lovld, hivld
         do j = lovld, hivld
            do k = lovld, hivld
               if (.not.
     &              (i .ge. d2lo .and. i .le. d2hi) .and.
     &              (j .ge. d2lo .and. j .le. d2hi) .and.
     &              (k .ge. d2lo .and. k .le. d2hi)) then
                  exhr = phiex(i, j, k)
                  calchr = phi(i, j, k)
                  dif = calchr - exhr
                  difabs = dabs(dif)
                  clcabs = dabs(calchr)
                  exabs = dabs(exhr)
                  sum1e = sum1e + exabs
                  sum1c = sum1c + clcabs
                  sum1d = sum1d + difabs
                  sum2e = sum2e + (exabs*exabs)
                  sum2c = sum2c + (clcabs*clcabs)
                  sum2d = sum2d + (difabs*difabs)
                  if (difabs .gt. difmax) then
                     difmax = difabs
                     difi = i
                     difj = j
                     difk = k
                  endif
                  if (exabs .gt. exmax) then
                     exmax = exabs
                     exi = i
                     exj = j
                     exk = k
                  endif
                  if (clcabs .gt. clcmax) then
                     clcmax = clcabs
                     clci = i
                     clcj = j
                     clck = k
                  endif
                  if (output) then
                     write (2,*) i, j, k, exhr, calchr, calchr-exhr
                  endif
               endif
            enddo
         enddo
      enddo
      dv = h**3
!     find 1-norms
      norm1e = sum1e * dv
      norm1c = sum1c * dv
      norm1d = sum1d * dv
!     find 2-norms
      norm2e = dsqrt(sum2e * dv)
      norm2c = dsqrt(sum2c * dv)
      norm2d = dsqrt(sum2d * dv)
      if (output) then
         close (unit=2)
      endif
      if (vrbos .gt. 0) then
         print *, 'NORMS of exact coarse solution, ',
     &        'calculated coarse solution, and difference'
         print *, 'displaying:  ',
     &        'point of maximum, max norm, L1 norm, L2 norm'
         print *, 'exact norms:  (', exi, exj, exk, ') ',
     &        exmax, norm1e, norm2e
         print *, 'calc norms:  (', clci, clcj, clck, ') ',
     &        clcmax, norm1c, norm2c
         print *, 'diff norms:  (', difi, difj, difk, ') ',
     &        difmax, norm1d, norm2d
      endif
      difnrm(0) = difmax
      difnrm(1) = norm1d
      difnrm(2) = norm2d

      return
      end

!=======================================================================

      subroutine precondrhs(rhs, l1, l2, l3, u1, u2, u3)
!     Replace rhs at interior points of its domain by
!     rhs + (h^2/12) * L_h(rhs)
!     where L_h is the 7-point discrete laplacian.
      implicit none
!     --- arguments ---
      integer l1, l2, l3, u1, u2, u3
      real*8 rhs(l1:u1, l2:u2, l3:u3)
!     --- locals ---
      integer i, j, k
      real*8 fac, rhscop(l1:u1, l2:u2, l3:u3)

!     add h^2/12 * L_7(rho)
!     = h^2/12 * (-6/h^2*rho(i) + sum_d (rho(i+e_d)/h^2 + rho(i-e_d)/h^2) )
!     = 1/12 * (-6 * rho(i) + sum_d (rho(i+e_d) + rho(i-e_d)) )
      fac = 1d0 / 12d0
      do k = l3+1, u3-1
         do j = l2+1, u2-1
            do i = l1+1, u1-1
               rhscop(i, j, k) = rhs(i, j, k) + fac *
     &              (-6d0 * rhs(i, j, k) +
     &              rhs(i+1, j, k) + rhs(i-1, j, k) +
     &              rhs(i, j+1, k) + rhs(i, j-1, k) +
     &              rhs(i, j, k+1) + rhs(i, j, k-1))
            enddo
         enddo
      enddo

      do k = l3+1, u3-1
         do j = l2+1, u2-1
            do i = l1+1, u1-1
               rhs(i, j, k) = rhscop(i, j, k)
            enddo
         enddo
      enddo

      return
      end
