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
      subroutine inhomodirichlet3d(f, lx, ly, lz, ux, uy, uz,
     &     dx, dy, dz, op, ier)
!     Nodal-point Poisson solver with inhomogeneous Dirichlet
!     boundary conditions on a cube.
!
!     n : input, number of unknowns in each direction, + 1.
!     n should be a highly composite number.  n must be even.
!     f(lx:ux, ly:uy, lz:uz): on input,
!     rhs in interior (lx+1:ux-1, ly+1:uy-1, lz+1:uz-1), and boundary
!     values on faces (points with at least one lower or upper coordinate);
!     on output, solution to Poisson's equation.
!     dx, dy, dz : input, mesh spacings.
!     op: input, operator type
!     (0=Spectral, 1=19-point Mehrstellen, 2=7point, 3=27-point Mehrstellen)
!     can't do op == 0
!     ier :  error code, 0 if and only if op code is legal
      implicit none
!     --- arguments ---
      integer lx, ly, lz, ux, uy, uz, op, ier
      real*8  f(lx:ux, ly:uy, lz:uz)
      real*8  dx, dy, dz
!     --- locals ---
      integer i, j, k, nx, ny, nz
      real*8 dx2i, dy2i, dz2i, c1, c2, c3, s2off, s3off
      real*8 h2inv, h(3)
      real*8 gxlo(lx:lx, ly:uy, lz:uz), gxhi(ux:ux, ly:uy, lz:uz)
      real*8 gylo(lx:ux, ly:ly, lz:uz), gyhi(lx:ux, uy:uy, lz:uz)
      real*8 gzlo(lx:ux, ly:uy, lz:lz), gzhi(lx:ux, ly:uy, uz:uz)

!     BACKGROUND
!     We want to solve:
!     L(phi) = rho
!     with boundary condition phi = phi_B on boundary.
!
!     Set phi = phi_I + phi_B
!     where phi_I is nonzero in interior, zero on boundary;
!     and phi_B is nonzero on boundary, zero in interior.
!     Since L is linear, we may rewrite:
!     L(phi_I) = rho - L(phi_B)
!     Then we have a homogeneous Dirichlet boundary condition for phi_I,
!     because phi_I was defined to be zero on the boundary.

!     In this subroutine, we save the initial boundary values in
!     gxlo, gxhi, gylo, gyhi, gzlo, gzhi;
!     then subtract L(phi_B) from rho (argument f).
!     Note that since phi_B = 0 in interior, L(phi_B) is nonzero only
!     at points adjacent to the boundary.

!     If the discrete operator is expressed as
!     L(f)[p] = sum_{d = 0:3} c_d * sum_{|q|=d} f(p + q)
!     then if e_b is unit vector pointing from p to the boundary face
!     and e1 and e2 are unit vectors in the other two coordinate
!     directions, then the only nonzero f's in the operator are
!     on the boundary, so this subroutine does:
!     f(p) := f(p) - (c_1 * f(p + e_b) + c_2 *
!     ( f(p+e_b+e1) + f(p+e_b-e1) + f(p+e_b+e2) + f(p+e_b-e2) ) + c_3 *
!     ( f(p+e_b-e1-e2) + f(p+e_b-e1+e2) + f(p+e_b+e1-e2) + f(p+e_b+e1+e2) ))

      dx2i = 1d0 / (dx * dx)
      dy2i = 1d0 / (dy * dy)
      dz2i = 1d0 / (dz * dz)
!     FIX this:  forget about h2inv and use each one separately
      if (dx .ne. dy .or. dx .ne. dz) then
         print *, 'cannot do anisotropic mesh spacings'
         call MAYDAY_ERROR()
      endif
      h2inv = dx2i

!      if (ux-lx .ne. uy-ly .or. ux-lx .ne. uz-lz) then
!         print *, 'cannot do anisotropic lengths'
!         call MAYDAY_ERROR()
!      endif

!     Save boundary values.
      do k = lz, uz
         do j = ly, uy
            gxlo(lx, j, k) = f(lx, j, k)
            gxhi(ux, j, k) = f(ux, j, k)
         enddo
      enddo
      do k = lz, uz
         do i = lx, ux
            gylo(i, ly, k) = f(i, ly, k)
            gyhi(i, uy, k) = f(i, uy, k)
         enddo
      enddo
      do j = ly, uy
         do i = lx, ux
            gzlo(i, j, lz) = f(i, j, lz)
            gzhi(i, j, uz) = f(i, j, uz)
         enddo
      enddo

      if (op .eq. 0) then
!     spectral discretization:  what is operator?
         print *, 'spectral operator not implemented'
         ier = 1
         return
      elseif (op .eq. 1) then
!     op 1 is 19-point operator
!     coefficients are:
!     -4 * h2inv for self
!     1/3 * h2inv for points differing by 1 in one coordinate
!     1/6 * h2inv for points differing by 1 in two coordinates
!     0 for points differing by 1 in three coordinates

!     Change f at points that are 1 away from faces.
!     There will be some repetition.
         c1 = h2inv / 3d0
         c2 = h2inv / 6d0
         do k = lz+1, uz-1
            do j = ly+1, uy-1
               s2off =
     &              f(lx, j-1, k  ) + f(lx, j+1, k  ) +
     &              f(lx, j  , k-1) + f(lx, j  , k+1)
               f(lx+1, j, k) = f(lx+1, j, k) -
     &              (c1*f(lx, j, k) + c2*s2off)

               s2off =
     &              f(ux, j-1, k  ) + f(ux, j+1, k  ) +
     &              f(ux, j  , k-1) + f(ux, j  , k+1)
               f(ux-1, j, k) = f(ux-1, j, k) -
     &              (c1*f(ux, j, k) + c2*s2off)
            enddo
         enddo
         do k = lz+1, uz-1
            do i = lx+1, ux-1
               s2off =
     &              f(i-1, ly, k  ) + f(i+1, ly, k  ) +
     &              f(i  , ly, k-1) + f(i  , ly, k+1)
               f(i, ly+1, k) = f(i, ly+1, k) -
     &              (c1*f(i, ly, k) + c2*s2off)

               s2off =
     &              f(i-1, uy, k  ) + f(i+1, uy, k  ) +
     &              f(i  , uy, k-1) + f(i  , uy, k+1)
               f(i, uy-1, k) = f(i, uy-1, k) -
     &              (c1*f(i, uy, k) + c2*s2off)

            enddo
         enddo
         do j = ly+1, uy-1
            do i = lx+1, ux-1
               s2off =
     &              f(i-1, j  , lz) + f(i+1, j  , lz) +
     &              f(i  , j-1, lz) + f(i  , j+1, lz)
               f(i, j, lz+1) = f(i, j, lz+1) -
     &              (c1*f(i, j, lz) + c2*s2off)

               s2off =
     &              f(i-1, j  , uz) + f(i+1, j  , uz) +
     &              f(i  , j-1, uz) + f(i  , j+1, uz)
               f(i, j, uz-1) = f(i, j, uz-1) -
     &              (c1*f(i, j, uz) + c2*s2off)
            enddo
         enddo
!     For places where you subtracted twice, add once again.
         do i = lx+1, ux-1
            f(i, ly+1, lz+1) = f(i, ly+1, lz+1) + c2*f(i, ly, lz)
            f(i, ly+1, uz-1) = f(i, ly+1, uz-1) + c2*f(i, ly, uz)
            f(i, uy-1, lz+1) = f(i, uy-1, lz+1) + c2*f(i, uy, lz)
            f(i, uy-1, uz-1) = f(i, uy-1, uz-1) + c2*f(i, uy, uz)
         enddo
         do j = ly+1, uy-1
            f(lx+1, j, lz+1) = f(lx+1, j, lz+1) + c2*f(lx, j, lz)
            f(lx+1, j, uz-1) = f(lx+1, j, uz-1) + c2*f(lx, j, uz)
            f(ux-1, j, lz+1) = f(ux-1, j, lz+1) + c2*f(ux, j, lz)
            f(ux-1, j, uz-1) = f(ux-1, j, uz-1) + c2*f(ux, j, uz)
         enddo
         do k = lz+1, uz-1
            f(lx+1, ly+1, k) = f(lx+1, ly+1, k) + c2*f(lx, ly, k)
            f(lx+1, uy-1, k) = f(lx+1, uy-1, k) + c2*f(lx, uy, k)
            f(ux-1, ly+1, k) = f(ux-1, ly+1, k) + c2*f(ux, ly, k)
            f(ux-1, uy-1, k) = f(ux-1, uy-1, k) + c2*f(ux, uy, k)
         enddo
      elseif (op .eq. 2) then
!     op 2 is 7-point operator (standard)
!     coefficients are:
!     -6 * h2inv for self
!     1 * h2inv for points differing by 1 in one coordinate
!     0 for points differing by 1 in two coordinates
!     0 for points differing by 1 in three coordinates

!     No points on intersections of faces will be used in operator.
!     Note that some points will be changed more than once, because
!     they are adjacent to more than one face.
         c1 = h2inv
         do k = lz+1, uz-1
            do j = ly+1, uy-1
               f(lx+1, j, k) = f(lx+1, j, k) - c1 * f(lx, j, k)
               f(ux-1, j, k) = f(ux-1, j, k) - c1 * f(ux, j, k)
            enddo
         enddo
         do k = lz+1, uz-1
            do i = lx+1, ux-1
               f(i, ly+1, k) = f(i, ly+1, k) - c1 * f(i, ly, k)
               f(i, uy-1, k) = f(i, uy-1, k) - c1 * f(i, uy, k)
            enddo
         enddo
         do j = ly+1, uy-1
            do i = lx+1, ux-1
               f(i, j, lz+1) = f(i, j, lz+1) - c1 * f(i, j, lz)
               f(i, j, uz-1) = f(i, j, uz-1) - c1 * f(i, j, uz)
            enddo
         enddo
      elseif (op .eq. 3) then
!     op 3 is 27-point operator
!     coefficients are:
!     -128/30 * h2inv for self
!     7/15 * h2inv for points differing by 1 in one coordinate
!     1/10 * h2inv for points differing by 1 in two coordinates
!     1/30 * h2inv for points differing by 1 in three coordinates

!     Change f at points that are 1 away from faces.
!     There will be some repetition.
         c1 = h2inv * 7d0/ 15d0
         c2 = h2inv / 10d0
         c3 = h2inv / 30d0
         do k = lz+1, uz-1
            do j = ly+1, uy-1
               s2off =
     &              f(lx, j-1, k  ) + f(lx, j+1, k  ) +
     &              f(lx, j  , k-1) + f(lx, j  , k+1)
               s3off =
     &              f(lx, j-1, k-1) + f(lx, j-1, k+1) +
     &              f(lx, j+1, k-1) + f(lx, j+1, k+1)
               f(lx+1, j, k) = f(lx+1, j, k) -
     &              (c1*f(lx, j, k) + c2*s2off + c3*s3off)

               s2off =
     &              f(ux, j-1, k  ) + f(ux, j+1, k  ) +
     &              f(ux, j  , k-1) + f(ux, j  , k+1)
               s3off =
     &              f(ux, j-1, k-1) + f(ux, j-1, k+1) +
     &              f(ux, j+1, k-1) + f(ux, j+1, k+1)
               f(ux-1, j, k) = f(ux-1, j, k) -
     &              (c1*f(ux, j, k) + c2*s2off + c3*s3off)
            enddo
         enddo
         do k = lz+1, uz-1
            do i = lx+1, ux-1
               s2off =
     &              f(i-1, ly, k  ) + f(i+1, ly, k  ) +
     &              f(i  , ly, k-1) + f(i  , ly, k+1)
               s3off =
     &              f(i-1, ly, k-1) + f(i-1, ly, k+1) +
     &              f(i+1, ly, k-1) + f(i+1, ly, k+1)
               f(i, ly+1, k) = f(i, ly+1, k) -
     &              (c1*f(i, ly, k) + c2*s2off + c3*s3off)

               s2off =
     &              f(i-1, uy, k  ) + f(i+1, uy, k  ) +
     &              f(i  , uy, k-1) + f(i  , uy, k+1)
               s3off =
     &              f(i-1, uy, k-1) + f(i-1, uy, k+1) +
     &              f(i+1, uy, k-1) + f(i+1, uy, k+1)
               f(i, uy-1, k) = f(i, uy-1, k) -
     &              (c1*f(i, uy, k) + c2*s2off + c3*s3off)
            enddo
         enddo
         do j = ly+1, uy-1
            do i = lx+1, ux-1
               s2off =
     &              f(i-1, j  , lz) + f(i+1, j  , lz) +
     &              f(i  , j-1, lz) + f(i  , j+1, lz)
               s3off =
     &              f(i-1, j-1, lz) + f(i-1, j+1, lz) +
     &              f(i+1, j-1, lz) + f(i+1, j+1, lz)
               f(i, j, lz+1) = f(i, j, lz+1) -
     &              (c1*f(i, j, lz) + c2*s2off + c3*s3off)

               s2off =
     &              f(i-1, j  , uz) + f(i+1, j  , uz) +
     &              f(i  , j-1, uz) + f(i  , j+1, uz)
               s3off =
     &              f(i-1, j-1, uz) + f(i-1, j+1, uz) +
     &              f(i+1, j-1, uz) + f(i+1, j+1, uz)
               f(i, j, uz-1) = f(i, j, uz-1) -
     &              (c1*f(i, j, uz) + c2*s2off + c3*s3off)
            enddo
         enddo
!     For places where you subtracted twice, add once again.
         do i = lx+1, ux-1
            f(i, ly+1, lz+1) = f(i, ly+1, lz+1) + c2*f(i, ly, lz)
     &           + c3*(f(i-1, ly, lz) + f(i+1, ly, lz))
            f(i, ly+1, uz-1) = f(i, ly+1, uz-1) + c2*f(i, ly, uz)
     &           + c3*(f(i-1, ly, uz) + f(i+1, ly, uz))
            f(i, uy-1, lz+1) = f(i, uy-1, lz+1) + c2*f(i, uy, lz)
     &           + c3*(f(i-1, uy, lz) + f(i+1, uy, lz))
            f(i, uy-1, uz-1) = f(i, uy-1, uz-1) + c2*f(i, uy, uz)
     &           + c3*(f(i-1, uy, uz) + f(i+1, uy, uz))
         enddo
         do j = ly+1, uy-1
            f(lx+1, j, lz+1) = f(lx+1, j, lz+1) + c2*f(lx, j, lz)
     &           + c3*(f(lx, j-1, lz) + f(lx, j+1, lz))
            f(lx+1, j, uz-1) = f(lx+1, j, uz-1) + c2*f(lx, j, uz)
     &           + c3*(f(lx, j-1, lz) + f(ux, j+1, lz))
            f(ux-1, j, lz+1) = f(ux-1, j, lz+1) + c2*f(ux, j, lz)
     &           + c3*(f(ux, j-1, lz) + f(ux, j+1, lz))
            f(ux-1, j, uz-1) = f(ux-1, j, uz-1) + c2*f(ux, j, uz)
     &           + c3*(f(ux, j-1, uz) + f(ux, j+1, uz))
         enddo
         do k = lz+1, uz-1
            f(lx+1, ly+1, k) = f(lx+1, ly+1, k) + c2*f(lx, ly, k)
     &           + c3*(f(lx, ly, k-1) + f(lx, ly, k+1))
            f(lx+1, uy-1, k) = f(lx+1, uy-1, k) + c2*f(lx, uy, k)
     &           + c3*(f(lx, uy, k-1) + f(lx, uy, k+1))
            f(ux-1, ly+1, k) = f(ux-1, ly+1, k) + c2*f(ux, ly, k)
     &           + c3*(f(ux, ly, k-1) + f(ux, ly, k+1))
            f(ux-1, uy-1, k) = f(ux-1, uy-1, k) + c2*f(ux, uy, k)
     &           + c3*(f(ux, uy, k-1) + f(ux, uy, k+1))
         enddo
!     corrections to the corrections
         f(lx+1, ly+1, lz+1) = f(lx+1, ly+1, lz+1) - c3 * f(lx, ly, lz)
         f(lx+1, ly+1, uz-1) = f(lx+1, ly+1, uz-1) - c3 * f(lx, ly, uz)
         f(lx+1, uy-1, lz+1) = f(lx+1, uy-1, lz+1) - c3 * f(lx, uy, lz)
         f(lx+1, uy-1, uz-1) = f(lx+1, uy-1, uz-1) - c3 * f(lx, uy, uz)
         f(ux-1, ly+1, lz+1) = f(ux-1, ly+1, lz+1) - c3 * f(ux, ly, lz)
         f(ux-1, ly+1, uz-1) = f(ux-1, ly+1, uz-1) - c3 * f(ux, ly, uz)
         f(ux-1, uy-1, lz+1) = f(ux-1, uy-1, lz+1) - c3 * f(ux, uy, lz)
         f(ux-1, uy-1, uz-1) = f(ux-1, uy-1, uz-1) - c3 * f(ux, uy, uz)
      else
         print *, 'inhomopoisson3d:  unknown op code ', op
         ier = 2
         return
      endif

!     Call Poisson solver with homogeneous boundary conditions.

      nx = ux - lx
      ny = uy - ly
      nz = uz - lz
      h(1) = dx
      h(2) = dy
      h(3) = dz
      call poisson3d(f, nx, ny, nz, h, op)

!     Fill in saved boundary conditions.

      do k = lz, uz
         do j = ly, uy
            f(lx, j, k) = gxlo(lx, j, k)
            f(ux, j, k) = gxhi(ux, j, k)
         enddo
      enddo
      do k = lz, uz
         do i = lx, ux
            f(i, ly, k) = gylo(i, ly, k)
            f(i, uy, k) = gyhi(i, uy, k)
         enddo
      enddo
      do j = ly, uy
         do i = lx, ux
            f(i, j, lz) = gzlo(i, j, lz)
            f(i, j, uz) = gzhi(i, j, uz)
         enddo
      enddo

      ier = 0

      return
      end
