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
      subroutine evaldirectfaces(ev, q,
     &     n1, s1, s2, n2c, layer, h, vrbos, ier)
      implicit none
!     --- arguments ---
      integer n1, s1, s2, n2c, layer
!     ev is at points of faces of D2, plus a layer
!     petermc, 8 Jan 2004:  last two dimensions are x/y/z and lo/hi
      real*8 ev(-layer:n2c+layer, -layer:n2c+layer, 1:3, 0:1)
!     q is at points of faces of D1
!     petermc, 8 Jan 2004:  last two dimensions are x/y/z and lo/hi
      real*8 q(0:n1, 0:n1, 1:3, 0:1)
      real*8 h
      integer vrbos
      integer ier
!     --- locals ---
      integer ncrsp, idims, idimd, isides, isided
      integer iev, icp, ierev

!     ncrsp is length of the first two dimensions of ev
      ncrsp = n2c + 2 * layer

!     loop over destination faces
      do isided = 0, 1
         do idimd = 1, 3
!     zero out the evaluations first
            call zerosquare(ev(-layer, -layer, idimd, isided), ncrsp)
            iev = (2*isided - 1) * idimd
!     loop over source faces
            do isides = 0, 1
               do idims = 1, 3
                  icp = (2*isides - 1) * idims
                  call evaldirect(ev(-layer, -layer, idimd, isided),
     &                 q(0, 0, idims, isided),
     &                 iev, icp,
     &                 n1, s1, s2, n2c, layer, h, vrbos, ierev)
                  if (ierev . ne. 0) then
                     ier = ierev
                     return
                  endif
               enddo
            enddo
!     end of loop over source faces
         enddo
      enddo
!     end of loop over destination faces

      ier = 0

      return
      end

!=======================================================================

      subroutine evaldirect(ev, q, iev, icp,
     &     n1, s1, s2, n2c, layer, h, vrbos, ier)
      implicit none
!     --- arguments ---
      integer iev, icp, n1, s1, s2, n2c, layer
!     ev is at coarse points of a face of D2, plus a layer
      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
!     q is at points of a face of D1
      real*8 q(0:n1, 0:n1)
      real*8 h
      integer vrbos
      integer ier
!     --- locals ---
      integer plsrc, sidsrc, isrc1, isrc2, src1, src2
      integer pldst, siddst, idst1, idst2, dest1, dest2
      integer n, n2, ptchsz
      integer ierg, iergi, intrul
      real*8 sum, xsrc(3), xdst(3), r, pi, hpatch
      real*8 bsrc(3), bdst(3), fac, newtrm, xdif(3), wt(0:n1)
      real*8 dsqrt, datan

      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, ierg)

      call getplaneindices(iev, pldst, siddst, idst1, idst2, ierg)

!     length of D
      n = n1 - 2*s1
!     length of D2
      n2 = n1 + 2*s2
!     coarsening factor
      ptchsz = n2 / n2c

!     Set bsrc(1:3) to be space coordinates of
!     the point indexed by (0, 0) on the source face of D1,
!     where indices are in the range (0:n1, 0:n1).
!     In the dimensions isrc1 and isrc2, a face of D1 extends
!     in space from -s1*h to (n+s1)*h = (n1-s1)*h.
!     Space coordinates of point indexed by (isrc1, isrc2) are
!     xsrc(isrc1) = bsrc(isrc1) + isrc1 * h
!     xsrc(isrc2) = bsrc(isrc2) + isrc2 * h
!     xsrc(plsrc) = bsrc(plsrc)
      bsrc(isrc1) = -s1 * h
      bsrc(isrc2) = -s1 * h
      if (sidsrc .eq. -1) then
         bsrc(plsrc) = -s1 * h
      elseif (sidsrc .eq. 1) then
         bsrc(plsrc) = (n + s1) * h
      else
         print *, 'illegal sidsrc'
      endif
      if (vrbos .gt. 0) then
         print *, 'source face ', icp,
     &        ' has (0,0) at ', bsrc(1), bsrc(2), bsrc(3),
     &        ', fixed dimension ', plsrc
      endif

!     Set bdst(1:3) to be space coordinates of
!     the point (0, 0) on the source face of D2,
!     where indices are in the range (-layer:n2c+layer, -layer:n2c+layer).
!     In the dimensions idst1 and idst2, a face of UNextended D2 extends
!     from (-s1-s2)*h to (n+s1+s2)*h = (n1-s1+s2)*h.
      bdst(idst1) = (-s1-s2)*h
      bdst(idst2) = (-s1-s2)*h
      if (siddst .eq. -1) then
         bdst(pldst) = (-s1 - s2) * h
      elseif (siddst .eq. 1) then
         bdst(pldst) = (n + s1 + s2) * h
      else
         print *, 'illegal siddst'
      endif
      if (vrbos .gt. 0) then
         print *, 'dest face ', iev,
     &        ' has (0,0) at ', bdst(1), bdst(2), bdst(3),
     &        ', fixed dimension ', pldst
      endif

      pi = 4 * datan(1d0)
      fac = 0.25d0 / pi

!     1 trapezoidal, 2 Simpson, 4 Bode, 8 Newton-Cotes-nine
      intrul = 4
!     Bode's rule here (intrul=4) is an improvement on
!     Simpson's (intrul=2) when n is big
!     -- can get 1/30 of minimum for n=512 and far faces.
!     Changing to intrul = 8 doesn't buy you anything.
      call getintegral1dweights(wt, 0, n1, h, intrul, iergi)

!     hpatch is spacing between coarse points on bound(D2).
      hpatch = ptchsz * h
!     xsrc(plsrc) and xdst(pldst) are fixed
      xsrc(plsrc) = bsrc(plsrc)
      xdst(pldst) = bdst(pldst)

!     loop over evaluation points, which are coarse points of bound(D2).
      do dest2 = -layer, n2c + layer
         xdst(idst2) = bdst(idst2) + dest2 * hpatch
         do dest1 = -layer, n2c + layer
            xdst(idst1) = bdst(idst1) + dest1 * hpatch
!     xdif(1:3) is indexed by source coordinates
            xdif(plsrc) = xdst(plsrc) - xsrc(plsrc)

!     Get old sum at destination point.
            sum = ev(dest1, dest2)
!     loop over source points
            do src2 = 0, n1
               xsrc(isrc2) = bsrc(isrc2) + src2 * h
               xdif(isrc2) = xdst(isrc2) - xsrc(isrc2)

               do src1 = 0, n1
                  xsrc(isrc1) = bsrc(isrc1) + src1 * h
                  xdif(isrc1) = xdst(isrc1) - xsrc(isrc1)

!                  r = dsqrt((xdst(1) - xsrc(1))**2 +
!     &                 (xdst(2) - xsrc(2))**2 +
!     &                 (xdst(3) - xsrc(3))**2)

                  r = dsqrt(xdif(1)*xdif(1) +
     &                 xdif(2)*xdif(2) +
     &                 xdif(3)*xdif(3))

                  newtrm = fac * wt(src1) * wt(src2) *
     &                 q(src1, src2) / r

!     ev(dest1, dest2) = ev(dest1, dest2) + newtrm
!     It's more efficient to calculate local sum and then later
!     add to ev(dest1, dest2).
                  sum = sum + newtrm
               enddo
            enddo
!     end of loop over source points

!     Recall that we had set sum = ev(dest1, dest2) before adding new terms.
            ev(dest1, dest2) = sum
!     print *, dest1, dest2, sum

         enddo
      enddo
!     end of loop over evaluation points

      ier = 0

      return
      end
