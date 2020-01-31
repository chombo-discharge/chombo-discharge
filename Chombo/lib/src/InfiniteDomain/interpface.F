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
      subroutine interpfaces(fine, coarse, intrpc,
     &     nfine, ptchsz, ncoar, layer)
!     Interpolate from coarse points to fine points on all faces of a cube.
!
!     fine : output, interpolated data on all faces,
!            with faces indexed by dimension 1:3 and low/high 0:1
!     coarse : input, data at coarse points
!     intrpc : input, interpolation polynomial coefficients from
!              call getpolyinterpc(intrpc, ptchsz, layer)
!     nfine : input, dimensions of fine in number of cells
!     ptchsz : input, refinement ratio between coarse and fine
!     ncoar : input, nfine / ptchsz
!     layer : input, width of layer of extra points of coarse
!
!     This subroutine performs two sweeps of interpolation, one
!     in each direction along the face.
!     One-dimensional interpolation at each fine point is from
!     2*layer+2 coarse points, with a polynomial of degree 2*layer+1.
!
!     The error in the approximation is O(H^(2*layer+2))
!     where H is the coarse mesh spacing, H = ptchsz * h
!     where h is the fine mesh spacing.
!     If H = sqrt(h) then error is O(h^(layer+1)).
      implicit none
!     --- arguments ---
      integer nfine, ptchsz, ncoar, layer
!     petermc, 8 Jan 2004:  last two dimensions are x/y/z and lo/hi
      real*8 coarse(-layer:ncoar+layer, -layer:ncoar+layer, 1:3, 0:1)
      real*8 intrpc(1:ptchsz-1, -layer:layer+1)
      real*8 fine(0:nfine, 0:nfine, 1:3, 0:1)
!     --- locals ---
      integer idim, iside

      do iside = 0, 1
         do idim = 1, 3
            call interpface(fine(0, 0, idim, iside),
     &           coarse(-layer, -layer, idim, iside), intrpc,
     &           nfine, nfine, ptchsz, ncoar, ncoar, layer)
         enddo
      enddo

      return
      end


      subroutine interpface(fine, coarse, intrpc,
     &     n1, n2, ptchsz, n1c, n2c, layer)
!     Interpolate from coarse points to fine points on a rectangular face.
!
!     fine : output, interpolated data on a face
!     coarse : input, data at coarse points
!     intrpc : input, interpolation polynomial coefficients from
!              call getpolyinterpc(intrpc, ptchsz, layer)
!     n1, n2 : input, dimensions of fine in number of cells
!     ptchsz : input, refinement ratio between coarse and fine
!     n1c : input, n1 / ptchsz
!     n2c : input, n2 / ptchsz
!     layer : input, width of layer of extra points of coarse
!
!     This subroutine performs two sweeps of interpolation, one
!     in each direction along the face.
!     One-dimensional interpolation at each fine point is from
!     2*layer+2 coarse points, with a polynomial of degree 2*layer+1.
!
!     The error in the approximation is O(H^(2*layer+2))
!     where H is the coarse mesh spacing, H = ptchsz * h
!     where h is the fine mesh spacing.
!     If H = sqrt(h) then error is O(h^(layer+1)).
      implicit none
!     --- arguments ---
      integer n1, n2, ptchsz, n1c, n2c, layer
      real*8 coarse(-layer:n1c+layer, -layer:n2c+layer)
      real*8 intrpc(1:ptchsz-1, -layer:layer+1)
      real*8 fine(0:n1, 0:n2)
!     --- locals ---
      integer icrse, jcrse, ifine, jfine, coff, foff, offd1, offd2
      real*8 sum
!     extension of fine with second coordinate coarse
      real*8 fineex(0:n1, -layer:n2c+layer)
      logical vrbose

      vrbose = .false.

      if (vrbose) then
         print *, 'Coarse evaluation on a face of D2; coarse range',
     &        -layer, ' to ', n2c+layer, n1c+layer
         do icrse = -layer, n1c+layer
            do jcrse = -layer, n2c+layer
               print *, icrse, jcrse, coarse(icrse, jcrse)
            enddo
         enddo
      endif

      offd1 = (n1c * ptchsz - n1) / 2
      offd2 = (n2c * ptchsz - n2) / 2
!     offdst == 0 iff ptchsz divides n2

!     Example with offdst=0 (and layer=1, ptchsz=4):
!     start with values at points marked 0,
!     then in first pass fill in those at points marked 1,
!     and in second pass fill in those at points marked 2.
!
!     0   0111011101110111011101110   0
!
!
!
!     0   0111011101110111011101110   0
!         2222222222222222222222222
!         2222222222222222222222222
!         2222222222222222222222222
!     0   0111011101110111011101110   0
!         2222222222222222222222222
!         2222222222222222222222222
!         2222222222222222222222222
!     0   0111011101110111011101110   0
!
!
!
!     0   0111011101110111011101110   0

!     first pass:  interpolate in first dimension
!     from points (ptchsz*(0:n2c), ptchsz*(-layer:n2c+layer))
!     to points (0:n2, ptchsz*(-layer:n2c+layer))

!     without a call to donothing, error goes way up
      call donothing

      do jcrse = -layer, n2c + layer
         do ifine = 0, n1
!     add ptchsz because if ifine-offdst < 0 then mod < 0.
!     note 0 <= offdst < ptchsz.
            foff = mod(ifine - offd1 + ptchsz, ptchsz)
            icrse = (ifine - foff) / ptchsz
            if (foff .eq. 0) then
!     projection from (coarse, coarse)
               fineex(ifine, jcrse) = coarse(icrse, jcrse)
            else
!     interpolation of data at (fine, coarse) points
!     from (coarse, coarse) points offset by multiples of (ptchsz, 0)
               sum = 0d0
               do coff = -layer, layer+1
                  sum = sum + intrpc(foff, coff) *
     &                 coarse(icrse + coff, jcrse)
               enddo
               fineex(ifine, jcrse) = sum
            endif
         enddo
      enddo

!     second pass:  interpolate in second dimension
!     from points (0:n2, ptchsz*(-layer:n2c+layer))
!     to points (0:n2, 0:n2)
      do jfine = 0, n2
         foff = mod(jfine - offd2 + ptchsz, ptchsz)
         jcrse = (jfine - foff) / ptchsz
         if (foff .eq. 0) then
            do ifine = 0, n1
!     projection from (fine, coarse)
               fine(ifine, jfine) = fineex(ifine, jcrse)
            enddo
         else
!     interpolation of data at (fine, fine) points
!     from (fine, coarse) points offset by multiples of (0, ptchsz)
            do ifine = 0, n1
               sum = 0d0
               do coff = -layer, layer+1
                  sum = sum + intrpc(foff, coff) *
     &                 fineex(ifine, jcrse + coff)
               enddo
               fine(ifine, jfine) = sum
            enddo
         endif
      enddo

      if (vrbose) then
         print *, 'Evaluation on a face.'
         do ifine = 0, n2
            do jfine = 0, n1
               print *, ifine, jfine, fine(ifine, jfine)
            enddo
         enddo
      endif

      return
      end


      subroutine donothing
      return
      end
