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

#include "REAL.H"
#include "CONSTANTS.H"

!=======================================================================

      integer function gcf2(n1, n2)
!     Returns greatest common factor of nonnegative integers n1 and n2.
      implicit none
!     --- arguments ---
      integer n1, n2
!     --- locals ---
      integer gcfnow, y, quot, r

      gcfnow = n1
      y = n2
      do while (y .gt. 0)
         quot = gcfnow / y
         r = gcfnow - y * quot
         gcfnow = y
         y = r
      enddo
      gcf2 = gcfnow

      return
      end

!=======================================================================

      integer function gcfarray(arr, n)
!     Returns greatest common factor of the n nonnegative integers
!     in array arr.
      implicit none
!     --- arguments ---
      integer n, arr(n)
!     --- locals ---
      integer gcfnow, i
      integer gcf2

      gcfnow = arr(1)
      do i = 2, n
         gcfnow = gcf2(gcfnow, arr(i))
      enddo

      gcfarray = gcfnow

      return
      end

!=======================================================================

      integer function minarray(arr, n)
!     Returns minimum of the n integers in array arr.
      implicit none
!     --- arguments ---
      integer n, arr(n)
!     --- locals ---
      integer minsf, i

      minsf = arr(1)
      do i = 2, n
         if (arr(i) .lt. minsf) then
            minsf = arr(i)
         endif
      enddo

      minarray = minsf

      return
      end

!=======================================================================

      integer function ceil(a)
      implicit none
!     --- arguments ---
      REAL_T a
!     --- locals ---
      integer i

      i = int(a)
      if (a .gt. i) then
         i = i + 1
      endif

      ceil = i

      return
      end

!=======================================================================

      subroutine fillextend(big, n, layer, small)
!     Copy a small 3-D array into a bigger one, and fill the bigger one
!     with zeroes outside the range of the smaller one.
!
!     big:  output, the big array, indices -layer:n+layer in each dimension
!     n:  input, length in cells of the smaller array
!     layer:  input, width of buffer between big and small arrays
!     small:  input, the small array, indices 0:n in each dimension
      implicit none
!     --- arguments ---
      integer n, layer
      REAL_T big(-layer:n+layer, -layer:n+layer, -layer:n+layer)
      REAL_T small(0:n, 0:n, 0:n)
!     --- locals ---
      integer i, j, k

      do k = -layer, n+layer
         do j = -layer, n+layer
            do i = -layer, n+layer
               big(i, j, k) = zero
            enddo
         enddo
      enddo
      do k = 0, n
         do j = 0, n
            do i = 0, n
               big(i, j, k) = small(i, j, k)
            enddo
         enddo
      enddo

      return
      end

!=======================================================================

      subroutine zerosquare(arr, n)
!     Fill a square 2-D array with zeroes.
!
!     arr:  output, the array, indices 0:n in both dimensions
!     n:  input, length of array
      implicit none
!     --- arguments ---
      integer n
      REAL_T arr(0:n, 0:n)
!     --- locals ---
      integer i, j

      do j = 0, n
         do i = 0, n
            arr(i, j) = zero
         enddo
      enddo

      return
      end

!=======================================================================

      subroutine zerosquarefaces(arr, n)
!     Fill an array representing face data with zeroes.
!
!     arr:  output, the array, indices 0:n in both dimensions on each face
!           and faces indexed by dimension 1:3 and low/high sides 0:1
!     n:  input, length of array
      implicit none
!     --- arguments ---
      integer n
      REAL_T arr(0:n, 0:n, 1:3, 0:1)
!     --- locals ---
      integer i, j, idim, iside

      do iside = 0, 1
         do idim = 1, 3
            do j = 0, n
               do i = 0, n
                  arr(i, j, idim, iside) = zero
               enddo
            enddo
         enddo
      enddo

      return
      end

!=======================================================================

      subroutine zerocube(arr, n)
!     Fill a cubic 3-D array with zeroes.
!
!     arr:  output, the array, indices 0:n in all dimensions
!     n:  input, length of array
      implicit none
!     --- arguments ---
      integer n
      REAL_T arr(0:n, 0:n, 0:n)
!     --- locals ---
      integer i, j, k

      do k = 0, n
         do j = 0, n
            do i = 0, n
               arr(i, j, k) = zero
            enddo
         enddo
      enddo

      return
      end

!=======================================================================

      subroutine getpolyinterpc(intrpc, ptchsz, layer)
!     Get coefficients of Lagrange interpolating polynomials
!     to interpolate values at fine points from coarse points.
!
!     intrpc:  output, the polynomial coefficients
!     ptchsz:  input, coarse/fine refinement ratio
!     layer:  input, width of the layer of coarse points beyond the
!             immediate neighbors to be used in interpolation
!
!     From values at integer points f(-layer), ..., f(layer+1)
!     we interpolate to the point i/ptchsz (1 <= i <= ptchsz-1) by:
!     f(i/ptchsz) = sum_{j=-layer:layer+1} intrpc(i, j) * f(j)
!
!     The approximation uses 2*layer+2 coarse points and hence a
!     polynomial of degree 2*layer+1 and so approximation error is
!     O(H^(2*layer+2)), where H is the coarse mesh spacing,
!     H = ptchsz * h, where h is the fine mesh spacing.
!     If H = sqrt(h) then approximation error is O(h^(layer+1)).
      implicit none
!     --- arguments ---
      integer ptchsz, layer
      REAL_T intrpc(1:ptchsz-1, -layer:layer+1)
!     --- locals ---
      integer ifine, icrse, icrso
      REAL_T xfine, xcrse, xcrso, prd

      do ifine = 1, ptchsz-1
         xfine = dfloat(ifine) / dfloat(ptchsz)
         do icrse = -layer, layer+1
            xcrse = dfloat(icrse)
!     find coefficient of f(xcrse) at xfine
            prd = 1d0
            do icrso = -layer, layer+1
               if (icrso .ne. icrse) then
                  xcrso = dfloat(icrso)
                  prd = prd * (xfine - xcrso) / (xcrse - xcrso)
               endif
            enddo
            intrpc(ifine, icrse) = prd
         enddo
      enddo

      return
      end

!=======================================================================

      subroutine setavgedgesall(n, g)
!     Replace values at edges of all boundary faces of cube by their averages.
!
!     n:  input, length in cells of each face
!     g:  on input, values on all faces;
!         on output, edge points are replaced by means.
!         faces indexed by dimension 1:3 and low/high sides 0:1
      implicit none
!     --- arguments ---
      integer n
      REAL_T g(0:n, 0:n, 1:3, 0:1)

      call setavgedges(n, g(0, 0, 1, 0), g(0, 0, 1, 1),
     &     g(0, 0, 2, 0), g(0, 0, 2, 1),
     &     g(0, 0, 3, 0), g(0, 0, 3, 1))

      return
      end

!=======================================================================

      subroutine setavgedges(n, gxlo, gxhi, gylo, gyhi, gzlo, gzhi)
!     Replace values at edges of boundary faces of a cube by their averages.
!
!     n:  input, length in cells of face
!     gxlo, gxhi, gylo, gyhi, gzlo, gzhi:  on input, values at faces;
!                on output, values at edge points are replaced by means of the
!                original values at all faces containing the point.
      implicit none
!     --- arguments ---
      integer n
!     gxlo at (0, 0:n, 0:n), gxhi at (n, 0:n, 0:n)
!     gylo at (0:n, 0, 0:n), gyhi at (0:n, n, 0:n)
!     gzlo at (0:n, 0:n, 0), gzhi at (0:n, 0:n, n)
      REAL_T gxlo(0:n, 0:n), gxhi(0:n, 0:n)
      REAL_T gylo(0:n, 0:n), gyhi(0:n, 0:n)
      REAL_T gzlo(0:n, 0:n), gzhi(0:n, 0:n)
!     --- locals ---
      integer k

!     8 corners:
!     to get (a, b, c), average gx(a, c), gy(a, c), gz(a, b)
!     where we pick gx to be gxlo if a=0, and gxhi if a=n;
!     similarly for gy and gz.

!     (0, 0, 0)
      call settomean3(gxlo(0, 0), gylo(0, 0), gzlo(0, 0))
!     (0, 0, n)
      call settomean3(gxlo(0, n), gylo(0, n), gzhi(0, 0))
!     (0, n, 0)
      call settomean3(gxlo(n, 0), gyhi(0, 0), gzlo(0, n))
!     (0, n, n)
      call settomean3(gxlo(n, n), gyhi(0, n), gzhi(0, n))
!     (n, 0, 0)
      call settomean3(gxhi(0, 0), gylo(n, 0), gzlo(n, 0))
!     (n, 0, n)
      call settomean3(gxhi(0, n), gylo(n, n), gzhi(n, 0))
!     (n, n, 0)
      call settomean3(gxhi(n, 0), gyhi(n, 0), gzlo(n, n))
!     (n, n, n)
      call settomean3(gxhi(n, n), gyhi(n, n), gzhi(n, n))

!     12 edges:
!     to get (a1, a2, a3), average gd1, gd2
!     where d1 and d2 are the UNchanged components,
!     and we pick gd1 to be gd1lo if a[d1]=0, gd1hi if a[d1]=n;
!     similarly for gd2.

!     (0:n, 0, 0)
      do k = 1, n-1
         call settomean2(gylo(k, 0), gzlo(k, 0))
      enddo
!     (0:n, 0, n)
      do k = 1, n-1
         call settomean2(gylo(k, n), gzhi(k, 0))
      enddo
!     (0:n, n, 0)
      do k = 1, n-1
         call settomean2(gyhi(k, 0), gzlo(k, n))
      enddo
!     (0:n, n, n)
      do k = 1, n-1
         call settomean2(gyhi(k, n), gzhi(k, n))
      enddo

!     (0, 0:n, 0)
      do k = 1, n-1
         call settomean2(gxlo(k, 0), gzlo(0, k))
      enddo
!     (0, 0:n, n)
      do k = 1, n-1
         call settomean2(gxlo(k, n), gzhi(0, k))
      enddo
!     (n, 0:n, 0)
      do k = 1, n-1
         call settomean2(gxhi(k, 0), gzlo(n, k))
      enddo
!     (n, 0:n, n)
      do k = 1, n-1
         call settomean2(gxhi(k, n), gzhi(n, k))
      enddo

!     (0, 0, 0:n)
      do k = 1, n-1
         call settomean2(gxlo(0, k), gylo(0, k))
      enddo
!     (0, n, 0:n)
      do k = 1, n-1
         call settomean2(gxlo(n, k), gyhi(0, k))
      enddo
!     (n, 0, 0:n)
      do k = 1, n-1
         call settomean2(gxhi(0, k), gylo(n, k))
      enddo
!     (n, n, 0:n)
      do k = 1, n-1
         call settomean2(gxhi(n, k), gyhi(n, k))
      enddo

      return
      end

!=======================================================================

      subroutine settomean3(a, b, c)
!     Replace a, b, c by their mean.
      implicit none
!     --- arguments ---
      REAL_T a, b, c
!     --- locals ---
      REAL_T mean

      mean = (a + b + c) * third
      a = mean
      b = mean
      c = mean

      return
      end

!=======================================================================

      subroutine settomean2(a, b)
!     Replace a, b by their mean.
      implicit none
!     --- arguments ---
      REAL_T a, b
!     --- locals ---
      REAL_T mean

      mean = (a + b) * half
      a = mean
      b = mean

      return
      end

!=======================================================

      REAL_T function evalpoly(coeffs, deg, x)
!     Returns the evaluation of the polynomial
!     sum_{d=0:deg} coeffs(d)*(x**d)
!     coeffs:  input array, coefficients of polynomial
!     deg:  input degree of polynomial
!     x:  input, point of evaluation of polynomial
      implicit none
!     --- arguments ---
      integer deg
      REAL_T coeffs(0:deg), x
!     --- locals ---
      integer d
      REAL_T sum

      sum = coeffs(deg)
      do d = deg-1, 0, -1
         sum = x * sum + coeffs(d)
      enddo
      evalpoly = sum

      return
      end

!=======================================================================

      subroutine getbuffersize(n, bufmin, s, ptchsz, ier)
!     Get size of buffer and size of patches to use with a particular
!     size of domain, when D1 is larger than D.
!
!     n:  input, length in cells of domain D of rhs;
!         n must be divisible by 4
!     bufmin:  input, minimum value of s
!     s:  output, width of buffer between D and D2
!     ptchsz:  output, length in cells of patches on domain of size n
!     ier:  output, error code, 0 if and only if n is in correct range
      implicit none
!     --- arguments ---
      integer n, bufmin, s, ptchsz, ier
!     --- locals ---
      integer n2c, n2cmin, ceil

      if (mod(n, 4) .ne. 0) then
         print *, 'getbuffersize:  n must be divisible by 4'
         ier = 4
         return
      elseif (n .lt. 4) then
         print *, 'getbuffersize:  n must be at least 4, not ', n
         ier = 5
         return
      endif

!     Want ptchsz ~ sqrt(n) and divisible by 4.
      ptchsz = four * ceil(sqrt(n * one) * fourth)

!     Need s2 >= ptchsz * sqrt(2).  Also s1 = s2.  And s = s1 + s2.
!     Hence s >= ptchsz * 2*sqrt(2).
!     D2 has length n + 2*(s1 + s2) = n + 2*s.
!     This length must be divisible by ptchsz.
!     Therefore n + 2*s = n2c * ptchsz.
!     And hence (n2c * ptchsz - n) / 2 = s >= ptchsz * 2*sqrt(2)
!     Therefore n2c >= n/ptchsz + 4*sqrt(2).
      n2c = ceil(4 * sqrt(two) + n/(ptchsz * one))

!     We have the condition s >= bufmin.
!     Use this to increase n2c if necessary.
      n2cmin = ceil((n + 2*bufmin) * one / (ptchsz * one))
      if (n2cmin .gt. n2c) then
         n2c = n2cmin
      endif

!     s = (n2c * ptchsz - n) / 2, and s must be even (because s1 = s2 = s/2).
!     But s is even because both ptchsz and n are divisible by 4.
      s = (n2c * ptchsz - n) / 2

      ier = 0

      return
      end

!=======================================================================

      subroutine get0buffersize(n, bufmin, s, ptchsz, ier)
!     Get size of buffer and size of patches to use with a particular
!     size of domain, when D1 == D.
!
!     n:  input, length in cells of domain D of rhs;
!         n must be divisible by 4
!     bufmin:  input, minimum value of s
!     s:  output, width of buffer between D=D1 and D2
!     ptchsz:  output, length in cells of patches on domain of size n
!     ier:  output, error code, 0 if and only if n is in correct range
      implicit none
!     --- arguments ---
      integer n, bufmin, s, ptchsz, ier
!     --- locals ---
      integer n2c, n2cmin, ceil

      if (mod(n, 4) .ne. 0) then
         print *, 'get0buffersize:  n must be divisible by 4'
         ier = 4
         return
      elseif (n .lt. 4) then
         print *, 'get0buffersize:  n must be at least 4, not ', n
         ier = 5
         return
      endif

!     Want ptchsz ~ sqrt(n) and divisible by 4.
      ptchsz = four * ceil(sqrt(n * one) * fourth)

!     Need s >= ptchsz * sqrt(2).
!     D2 has length n + 2*s.
!     This length must be divisible by ptchsz.
!     Therefore n + 2*s = n2c * ptchsz.
!     And hence (n2c * ptchsz - n) / 2 = s >= ptchsz * 2*sqrt(2)
!     Therefore n2c >= n/ptchsz + 4*sqrt(2).
      n2c = ceil(two * sqrt(two) + n/(ptchsz * one))

!     We have the condition s >= bufmin.
!     Use this to increase n2c if necessary.
      n2cmin = ceil((n + 2*bufmin) * one / (ptchsz * one))
      if (n2cmin .gt. n2c) then
         n2c = n2cmin
      endif

!     s = (n2c * ptchsz - n) / 2, and s must be even (because s1 = s2 = s/2).
!     But s is even because both ptchsz and n are divisible by 4.
      s = (n2c * ptchsz - n) / 2

      ier = 0

      return
      end

!=======================================================================

      subroutine getcbuffersize(n, bufmin, s, ptchsz, ier)
!     Get size of buffer and size of patches to use with a particular
!     size of domain, when D1 == D and ratio of buffer to n is constant.
!
!     n:  input, length in cells of domain D of rhs; n must be even and
!         in the range 16 to 2048
!     bufmin:  input, minimum value of s
!     s:  output, width of buffer between D=D1 and D2
!     ptchsz:  output, length in cells of patches on domain of size n
!     ier:  output, error code, 0 if and only if n is in correct range
      implicit none
!     --- arguments ---
      integer n, bufmin, s, ptchsz, ier
!     --- locals ---
      integer n2c, n2cmin, ceil
      REAL_T con

      if (mod(n, 4) .ne. 0) then
         print *, 'getcbuffersize:  n must be divisible by 4'
         ier = 4
         return
      elseif (n .lt. 4) then
         print *, 'getcbuffersize:  n must be at least 4, not ', n
         ier = 5
         return
      endif

      con = three / eight

!     Want ptchsz ~ sqrt(n) and divisible by 4.
      ptchsz = four * ceil(sqrt(n * one) * fourth)

      n2c = ((one + two * con) * n) / (ptchsz * one)

!     We have the condition s >= bufmin.
!     Use this to increase n2c if necessary.
      n2cmin = ceil((n + 2*bufmin) * one / (ptchsz * one))
      if (n2cmin .gt. n2c) then
         n2c = n2cmin
      endif

!     s = (n2c * ptchsz - n) / 2, and s must be even (because s1 = s2 = s/2).
!     But s is even because both ptchsz and n are divisible by 4.
      s = (n2c * ptchsz - n) / 2

      ier = 0

      return
      end


!=======================================================================

      subroutine get0buffersize3d(n3d, bufmin, s, ptchsz, ier)
!     Get size of buffer and size of patches to use with a particular
!     size of domain, when D1 == D.
!
!     n3d:  input, length in cells of sides of domain D of rhs;
!           each component must be divisible by 4
!     bufmin:  input, minimum value of s
!     s:  output, width of buffer between D=D1 and D2
!     ptchsz:  output, length in cells of patches on domain of size n
!     ier:  output, error code, 0 if and only if n is in correct range
      implicit none
!     --- arguments ---
      integer n3d(3), bufmin, s, ptchsz, ier
!     --- locals ---
      integer n2c(3), n2cmin, ceil, idim, n, pchsza(3), sall(3)
      integer ndif(3), gcfnd, psmin
      integer gcf2, gcfarray, minarray

      do idim = 1, 3
         n = n3d(idim)
         if (mod(n, 4) .ne. 0) then
            print *,
     &           'get0buffersize3d:  lengths must be divisible by 4'
            ier = 4
            return
         elseif (n .lt. 4) then
            print *,
     &           'get0buffersize3d:  lengths must be at least 4'
            ier = 5
            return
         endif
      enddo

!     Want ptchsz ~ sqrt(n) and divisible by 4.
      do idim = 1, 3
         n = n3d(idim)
         pchsza(idim) = four * ceil(sqrt(n * one) * fourth)
      enddo
!     We want patch size to divide all differences in n
      ndif(1) = abs(n3d(1) - n3d(2))
      ndif(2) = abs(n3d(1) - n3d(3))
      ndif(3) = abs(n3d(2) - n3d(3))
      gcfnd = gcfarray(ndif, 3)
      psmin = minarray(pchsza, 3)
      ptchsz = gcf2(psmin, gcfnd)

!     Need s >= ptchsz * sqrt(2).
!     D2 has length n + 2*s.
!     This length must be divisible by ptchsz.
!     Therefore n + 2*s = n2c * ptchsz.
!     And hence (n2c * ptchsz - n) / 2 = s >= ptchsz * 2*sqrt(2)
!     Therefore n2c >= n/ptchsz + 4*sqrt(2).
      do idim = 1, 3
         n = n3d(idim)
         n2c(idim) = ceil(two * sqrt(two) + n/(ptchsz * one))

!     We have the condition s >= bufmin.
!     Use this to increase n2c if necessary.
         n2cmin = ceil((n + 2*bufmin) * one / (ptchsz * one))
         if (n2cmin .gt. n2c(idim)) then
            n2c(idim) = n2cmin
         endif

!     s = (n2c * ptchsz - n) / 2, and s must be even (because s1 = s2 = s/2).
!     But s is even because both ptchsz and n are divisible by 4.
         sall(idim) = (n2c(idim) * ptchsz - n) / 2
      enddo

      s = 0
      do idim = 1, 3
         if (sall(idim) .gt. s) then
            s = sall(idim)
         endif
      enddo

      ier = 0

      return
      end


!=======================================================================

      subroutine getbuffersize3d(n3d, bufmin, s, ptchsz, ier)
!     Get size of buffer and size of patches to use with a particular
!     size of domain, when D1 is larger than D.
!
!     n3d:  input, length in cells of sides of domain D of rhs;
!           each component must be divisible by 4
!     bufmin:  input, minimum value of s
!     s:  output, width of buffer between D=D1 and D2
!     ptchsz:  output, length in cells of patches on domain of size n
!     ier:  output, error code, 0 if and only if n is in correct range
      implicit none
!     --- arguments ---
      integer n3d(3), bufmin, s, ptchsz, ier
!     --- locals ---
      integer n2c(3), n2cmin, ceil, idim, n, pchsza(3), sall(3)
      integer ndif(3), gcfnd, psmin
      integer gcf2, gcfarray, minarray

      do idim = 1, 3
         n = n3d(idim)
         if (mod(n, 4) .ne. 0) then
            print *,
     &           'getbuffersize3d:  lengths must be divisible by 4'
            ier = 4
            return
         elseif (n .lt. 4) then
            print *,
     &           'getbuffersize3d:  lengths must be at least 4'
            ier = 5
            return
         endif
      enddo

!     Want ptchsz ~ sqrt(n) and divisible by 4.
      do idim = 1, 3
         n = n3d(idim)
         pchsza(idim) = four * ceil(sqrt(n * one) * fourth)
      enddo
!     We want patch size to divide all differences in n
      ndif(1) = abs(n3d(1) - n3d(2))
      ndif(2) = abs(n3d(1) - n3d(3))
      ndif(3) = abs(n3d(2) - n3d(3))
      gcfnd = gcfarray(ndif, 3)
      psmin = minarray(pchsza, 3)
      ptchsz = gcf2(psmin, gcfnd)

!     Need s2 >= ptchsz * sqrt(2).  Also s1 = s2.  And s = s1 + s2.
!     Hence s >= ptchsz * 2*sqrt(2).
!     D2 has length n + 2*(s1 + s2) = n + 2*s.
!     This length must be divisible by ptchsz.
!     Therefore n + 2*s = n2c * ptchsz.
!     And hence (n2c * ptchsz - n) / 2 = s >= ptchsz * 2*sqrt(2)
!     Therefore n2c >= n/ptchsz + 4*sqrt(2).
      do idim = 1, 3
         n = n3d(idim)
         n2c(idim) = ceil(4 * sqrt(two) + n/(ptchsz * one))

!     We have the condition s >= bufmin.
!     Use this to increase n2c if necessary.
         n2cmin = ceil((n + 2*bufmin) * one / (ptchsz * one))
         if (n2cmin .gt. n2c(idim)) then
            n2c(idim) = n2cmin
         endif

!     s = (n2c * ptchsz - n) / 2, and s must be even (because s1 = s2 = s/2).
!     But s is even because both ptchsz and n are divisible by 4.
         sall(idim) = (n2c(idim) * ptchsz - n) / 2
      enddo

      s = 0
      do idim = 1, 3
         if (sall(idim) .gt. s) then
            s = sall(idim)
         endif
      enddo

      ier = 0

      return
      end
