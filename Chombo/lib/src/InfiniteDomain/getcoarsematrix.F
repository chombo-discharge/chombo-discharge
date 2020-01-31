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

      subroutine getallcoarsematrix(sdcrse, ncpts,
     &     n1, s2, ptchsz, npchs, deg, h, refc, cadd,
     &     cpind, cfac, cpow, nterms, vrbos, ier)
      implicit none
!     --- arguments ---
      integer n1, s2, ptchsz, npchs, deg, cadd, refc
!     number of terms in expansions of all derivatives of 1/r
      integer nterms, ncpts
      real*8 h
!     indices, powers, and constant factors in terms in
!     expansions of derivatives of 1/r
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      integer vrbos, ier
      real*8 sdcrse(((deg+1)*(deg+2))/2, ncpts,
     &     0:npchs-1, 0:npchs-1, 6)
!     --- locals ---
      integer ierev
      integer ori, isides, idims, icp

      ori = 0
!     loop over source faces
      do isides = 0, 1
         do idims = 1, 3
            icp = (2*isides - 1) * idims

            ori = ori + 1

            call getcoarsematrix(sdcrse(1, 1, 0, 0, ori),
     &           icp,
     &           n1, s2, ptchsz, npchs, deg, h, refc, cadd,
     &           cpind, cfac, cpow, nterms, ncpts, vrbos, ierev)

            if (ierev .ne. 0) then
               print *, 'error in getallcoarsematrix:  ', ierev
               ier = ierev
               return
            endif

         enddo
      enddo
!     end of loop over source faces

      ier = 0

      return
      end

c=======================================================================
      subroutine getcoarsematrix(sdcrse, icp,
     &     n1, s2, ptchsz, npchs, deg, h, refc, cadd,
     &     cpind, cfac, cpow, nterms, ncpts, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     of D1 on the coarse points in a region outside of D2.
!
!     => inputs, <= outputs
!     sdcrse <= the matrix
!     icp => +/- 1,2,3; source face with multipole expansions
!     n1 => length of domain D1
!     s2 => offset from D1 to D2
!     ptchsz => length of patches
!     npchs => number of patches in each dimension
!     deg => degree of multipole expansions
!     h => mesh spacing
!     refc => refinement ratio between fine and coarse
!     cadd => length of extra layer of coarse points beyond D2
!     cpind, cfac, cpow => from subroutine deriv1overr, describing
!                         expressions for partial derivatives of 1/r
!     nterms => number of terms in cpind, cfac, cpow
!     ncpts => number of coarse points
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer icp, n1, s2, ptchsz, npchs, deg, refc, cadd
      real*8 h
      integer nterms, ncpts
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      real*8 sdcrse(((deg+1)*(deg+2))/2, ncpts,
     &     0:npchs-1, 0:npchs-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, i, j, k, nc
      integer n2, n2c
      integer plsrc, sidsrc
      integer isrc1, isrc2, term, subtrm
      integer d2lo, d2hi
      real*8 hbig, hpatch
      real*8 dst(3), offset
      real*8 xisrc1, xisrc2, xplsrc, r, r2inv, sum
      integer iergpi, dd, pt
      integer endsz
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)
      real*8 srcpl, pctr1(0:npchs-1), pctr2(0:npchs-1)

!     number of coefficients
      nc = ((deg+1)*(deg+2))/2
!     length of D2
      n2 = n1 + 2*s2
!     number of coarse points within D2 in each dimension
      n2c = n2 / refc
!     physical width of each patch
      hpatch = ptchsz * h
!     mesh spacing for coarsened level
      hbig = refc * h

!     indlo(term) = min {i: cpind(i) == term}
!     indhi(term) = max {i: cpind(i) == term}
      do term = 1, nterms
         indhi(cpind(term)) = term
      enddo
      do term = nterms, 1, -1
         indlo(cpind(term)) = term
      enddo

!     set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
!     the two perpendicular directions on the source face
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getcoarsematrix:  illegal icp = ', icp
         ier = 1
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc
 110     format (1x, 'source face ', a2, i2)
      endif

!     endsz is the length of the first or last patch.
      if (npchs .gt. 1) then
         endsz = ptchsz - (npchs * ptchsz - n1) / 2
!     Thus 2 * endsz + (npchs - 2) * ptchsz = n1
      else
         endsz = n1
      endif

!
!     In this illustration, npchs = 3.
!     Patch centers are indicated by X.
!
!            :<-- ptchsz  -->:               :<-- ptchsz  -->:
!                :<- endsz ->:<-- ptchsz  -->:<- endsz ->:   :
!            +---+-----------+---------------+-----------+---+
!            |   :           :               :           :   |
!      - - - + - +-----------+---------------+-----------+ - | - - -
!        ^   |   |           |               |           |   |   ^
!        |   |   |           |               |           |   |   |
!      endsz |   |     X     |       X       |     X     |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   |
!      - - - + - +-----------+---------------+-----------+ - +   |
!        ^   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!     ptchsz |   |     X     |       X       |     X     |   |   n1
!        |   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   |
!      - - - + - +-----------+---------------+-----------+ - +   |
!        ^   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!      endsz |   |     X     |       X       |     X     |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   v
!      - - - + - +-----------+---------------+-----------+ - + - - -
!            |   :           :               :           :   |
!            +---+-----------+---------------+-----------+---+
!                :                                       :
!                :<--------------- n1 ------------------>:
!                :                                       :

!     A face of D1 extends in space
!     from s2*h to (n1+s2)*h = (n2-s2)*h in the isrc1 and isrc2 directions,
!     and fixed at either s2*h or (n1+s2)*h = (n2-s2)*h in plsrc direction.

!     The center of the (0, 0) patch is set off from the bottom left corner by
!     endsz*h/2 in the isrc1 and isrc2 directions.

!     srcpl is fixed physical coordinate in direction normal to patch.
      if (sidsrc .eq. -1) then
         srcpl = s2 * h
      elseif (sidsrc .eq. 1) then
         srcpl = (n1 + s2) * h
      endif

!     pctr1(patch1) is physical coordinate of center of source patch
!     indexed by patch1 in first parallel direction.
!     Bottom of patch 0 is at s2 * h.
      call getcenters(pctr1, npchs, s2 * h, ptchsz, endsz, h)

!     pctr2(patch2) is physical coordinate of center of source patch
!     indexed by patch2 in second parallel direction.
!     Bottom of patch 0 is at s2 * h.
      call getcenters(pctr2, npchs, s2 * h, ptchsz, endsz, h)

!     D2 extends from 0 to n2*h in each dimension.
!     refc need not divide n2:  if it doesn't, then the coarse outer points
!     are outside the faces of D2.
!     Set offset to the distance in fine cells from the boundary of D2
!     to the first coarse point in D2.
      offset = (mod(n2, refc) / 2) * h

!     d2lo is the lowest coarse coordinate of a coarse point inside D2
!      d2lo = -offc
      d2lo = 0
!     d2hi is the highest coarse coordinate of a coarse point inside D2
!      d2hi = nbigc + offc
      d2hi = n2c

      pt = 0
!     Loop over all destination points
      do k = -cadd, n2c + cadd
         dst(3) = k * hbig + offset
         do j = -cadd, n2c + cadd
            dst(2) = j * hbig + offset
            do i = -cadd, n2c + cadd
!     Check that we are outside D2:
!     at least one coordinate must be outside the range d2lo:d2hi.
               if (.not. (
     &              (i .ge. d2lo .and. i .le. d2hi) .and.
     &              (j .ge. d2lo .and. j .le. d2hi) .and.
     &              (k .ge. d2lo .and. k .le. d2hi))) then

                  pt = pt + 1
!                  print *, 'pt', pt, ' is ', i, j, k

                  dst(1) = i * hbig + offset

!     Loop over all patches on source face.
                  do patch2 = 0, npchs-1
                     do patch1 = 0, npchs-1
!     Destination point is dst(1:3), and we add to ev(i, j, k).
!     x is displacement vector from src to dst.
!     The multipole expansions are based on powers of x(isrc1) and x(isrc2).
                        xisrc1 = dst(isrc1) - pctr1(patch1)
                        xisrc2 = dst(isrc2) - pctr2(patch2)
                        xplsrc = dst(plsrc) - srcpl
!     r is length of x.
                        r = dsqrt(xisrc1**2 + xisrc2**2 + xplsrc**2)

                        r2inv = 1d0 / (r*r)

!     pow1(k) == x(isrc1)**k, for k = 0:deg
!     pow2(k) == x(isrc2)**k, for k = 0:deg
!     powr(k) == 1/r**(2*k+1), for k = 0:deg
                        pow1(0) = 1d0
                        pow2(0) = 1d0
                        powr(0) = 1d0 / r
                        do dd = 1, deg
                           pow1(dd) = xisrc1 * pow1(dd-1)
                           pow2(dd) = xisrc2 * pow2(dd-1)
                           powr(dd) = powr(dd-1) * r2inv
                        enddo

                        do term = 1, nc
                           sum = 0d0
                           do subtrm = indlo(term), indhi(term)
                              sum = sum +
     &                             cfac(subtrm) *
     &                             pow1(cpow(1, subtrm)) *
     &                             pow2(cpow(2, subtrm)) *
     &                             powr(cpow(3, subtrm))
                           enddo

!                           srcdst(term, i, j, k, patch1, patch2) =
!     &                          sum
!     petermc, 14 oct 2004:  changed sum to -sum
                           sdcrse(term, pt, patch1, patch2) = -sum

                        enddo

!     Recall that we set sum = ev(dest1, dest2) before adding new terms.
!     ev(dest1, dest2) = sum
                     enddo
                  enddo
!     end of loop over patches
               endif
!     end of check whether in range d2lo:d2hi
            enddo
         enddo
      enddo
!     end of loop over destination points

      if (pt .ne. ncpts) then
         print *, 'getcoarsematrix: ran through ', pt,
     &        ' coarse points instead of ', ncpts
         ier = 2
         return
      endif

      ier = 0

      return
      end


c=======================================================================
      subroutine getrectcoarsematrix(sdcrse, icp, nsrc, ndst, refc,
     &     ptchsz, npchs1, npchs2,
     &     deg, h, cadd,
     &     cpind, cfac, cpow, nterms, ncpts, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     of D1 on the coarse points in a region outside of D2.
!
!     => inputs, <= outputs
!     sdcrse <= the matrix
!     icp => +/- 1,2,3; source face with multipole expansions
!     nsrc => vector, lengths in cells of domain D1
!     ndst => vector, lengths in cells of domain D2
!     refc => refinement ratio between fine and coarse
!     ptchsz => length of patches
!     npchs1 => number of patches in 1st direction on face of D1
!     npchs2 => number of patches in 2nd direction on face of D1
!     deg => degree of multipole expansions
!     h => vector, mesh spacing in three dimensions
!     cadd => length of extra layer of coarse points beyond D2
!     cpind, cfac, cpow => from subroutine deriv1overr, describing
!                         expressions for partial derivatives of 1/r
!     nterms => number of terms in cpind, cfac, cpow
!     ncpts => number of coarse points
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer icp, nsrc(3), ndst(3), refc
      integer ptchsz, npchs1, npchs2
      integer deg, cadd
      real*8 h(3)
      integer nterms, ncpts
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      real*8 sdcrse(((deg+1)*(deg+2))/2, ncpts,
     &     0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer idir, buf(3)
      integer patch1, patch2, i, j, k, nc
      integer ndstc(3)
      integer plsrc, sidsrc
      integer isrc1, isrc2, term, subtrm
      integer d2lo(3), d2hi(3)
      real*8 hbig(3)
      real*8 dst(3), offset(3)
      real*8 xisrc1, xisrc2, xplsrc, r, r2inv, sum
      integer iergpi, dd, pt
      integer endp1, endp2
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)
      real*8 srcpl, pctr1(0:npchs1-1), pctr2(0:npchs2-1)

!     number of coefficients
      nc = ((deg+1)*(deg+2))/2

!     indlo(term) = min {i: cpind(i) == term}
!     indhi(term) = max {i: cpind(i) == term}
      do term = 1, nterms
         indhi(cpind(term)) = term
      enddo
      do term = nterms, 1, -1
         indlo(cpind(term)) = term
      enddo

!     set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
!     the two perpendicular directions on the source face
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getrectcoarsematrix:  illegal icp = ', icp
         ier = 1
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc
 110     format (1x, 'source face ', a2, i2)
      endif

      do idir = 1, 3
         buf(idir) = (ndst(idir) - nsrc(idir)) / 2
      enddo

!     endp1 is the length of the first or last patch in the first direction.
!     endp2 is the length of the first or last patch in the second direction.
!      endp1 = nsrc(isrc1) - (npchs1 - 1) * ptchsz
!      endp2 = nsrc(isrc2) - (npchs2 - 1) * ptchsz
      if (npchs1 .gt. 1) then
         endp1 = ptchsz - (npchs1 * ptchsz - nsrc(isrc1)) / 2
!     Thus 2 * endp1 + (npchs1 - 2) * ptchsz = nsrc(isrc1)
      else
         endp1 = nsrc(isrc1)
      endif
      if (npchs2 .gt. 1) then
         endp2 = ptchsz - (npchs2 * ptchsz - nsrc(isrc2)) / 2
!     Thus 2 * endp2 + (npchs2 - 2) * ptchsz = nsrc(isrc2).
      else
         endp2 = nsrc(isrc2)
      endif

!
!     In this illustration, npchs = 3.
!     Patch centers are indicated by X.
!
!            :<-- ptchsz  -->:               :<-- ptchsz  -->:
!                :<- endsz ->:<-- ptchsz  -->:<- endsz ->:   :
!            +---+-----------+---------------+-----------+---+
!            |   :           :               :           :   |
!      - - - + - +-----------+---------------+-----------+ - | - - -
!        ^   |   |           |               |           |   |   ^
!        |   |   |           |               |           |   |   |
!      endsz |   |     X     |       X       |     X     |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   |
!      - - - + - +-----------+---------------+-----------+ - +   |
!        ^   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!     ptchsz |   |     X     |       X       |     X     |   |   n1
!        |   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   |
!      - - - + - +-----------+---------------+-----------+ - +   |
!        ^   |   |           |               |           |   |   |
!        |   |   |           |               |           |   |   |
!      endsz |   |     X     |       X       |     X     |   |   |
!        |   |   |           |               |           |   |   |
!        v   |   |           |               |           |   |   v
!      - - - + - +-----------+---------------+-----------+ - + - - -
!            |   :           :               :           :   |
!            +---+-----------+---------------+-----------+---+
!                :                                       :
!                :<--------------- n1 ------------------>:
!                :                                       :

!     A face of D1 extends in space
!     from s2*h to (n1+s2)*h = (n2-s2)*h in the isrc1 and isrc2 directions,
!     and fixed at either s2*h or (n1+s2)*h = (n2-s2)*h in plsrc direction.

!     The center of the (0, 0) patch is set off from the bottom left corner by
!     endsz*h/2 in the isrc1 and isrc2 directions.

!     srcpl is fixed physical coordinate in direction normal to patch.
      if (sidsrc .eq. -1) then
         srcpl = buf(plsrc) * h(plsrc)
      elseif (sidsrc .eq. 1) then
         srcpl = (nsrc(plsrc) + buf(plsrc)) * h(plsrc)
      endif

!     (pctr1(patch1), pctr2(patch2)) are physical coordinates of
!     source patch indexed by (patch1, patch2).
!     Bottom of patch (0, 0) is at
!     (buf(isrc1) * h(isrc1), buf(isrc2) * h(isrc2)).
      call getcenters(pctr1, npchs1, buf(isrc1) * h(isrc1),
     &     ptchsz, endp1, h(isrc1))
      call getcenters(pctr2, npchs2, buf(isrc2) * h(isrc2),
     &     ptchsz, endp2, h(isrc2))

!     D2 extends in space from 0 to ndst*h.
!     refc need not divide ndst:  if it doesn't, then the coarse outer points
!     are outside the faces of D2.
      do idir = 1, 3
!     mesh spacing for coarsened level
         hbig(idir) = refc * h(idir)
!     Set offset to the distance in fine cells from the boundary of D2
!     to the first coarse point in D2.
         offset(idir) = mod(ndst(idir) / 2, refc) * h(idir)
!     d2lo is the lowest coarse coordinate of a coarse point inside D2
         d2lo(idir) = 0
!     The numerator should be divisible by the denominator
!     because ndst(idir) should be even
         ndstc(idir) = (ndst(idir) - 2*mod(ndst(idir)/2, refc)) / refc
!     d2hi is the highest coarse coordinate of a coarse point inside D2.
         d2hi(idir) = ndstc(idir)
!         offlo(idir) = (refc * d2hi(idir) - ndst(idir)) / 2
!     Set bdst to be physical coordinates of
!     the point with evaluation at ev(0, 0, 0).
!     This point has coordinates -offlo.
!         bdst(idir) = -offlo(idir) * h(idir)
      enddo

      pt = 0
!     Loop over all destination points
      do k = -cadd, ndstc(3) + cadd
         dst(3) = k * hbig(3) + offset(3)
         do j = -cadd, ndstc(2) + cadd
            dst(2) = j * hbig(2) + offset(2)
            do i = -cadd, ndstc(1) + cadd
!     Check that we are outside D2:
!     at least one coordinate must be outside the range d2lo:d2hi.
               if (.not. (
     &              (i .ge. d2lo(1) .and. i .le. d2hi(1)) .and.
     &              (j .ge. d2lo(2) .and. j .le. d2hi(2)) .and.
     &              (k .ge. d2lo(3) .and. k .le. d2hi(3)))) then
                  pt = pt + 1

                  dst(1) = i * hbig(1) + offset(1)

!     Loop over all patches on source face.
                  do patch2 = 0, npchs2-1
                     do patch1 = 0, npchs1-1
!     Destination point is dst(1:3), and we add to ev(i, j, k).
!     x is displacement vector from src to dst.
!     The multipole expansions are based on powers of x(isrc1) and x(isrc2).
                        xisrc1 = dst(isrc1) - pctr1(patch1)
                        xisrc2 = dst(isrc2) - pctr2(patch2)
                        xplsrc = dst(plsrc) - srcpl
!     r is length of x.
                        r = dsqrt(xisrc1**2 + xisrc2**2 + xplsrc**2)

                        r2inv = 1d0 / (r*r)

!     pow1(k) == x(isrc1)**k, for k = 0:deg
!     pow2(k) == x(isrc2)**k, for k = 0:deg
!     powr(k) == 1/r**(2*k+1), for k = 0:deg
                        pow1(0) = 1d0
                        pow2(0) = 1d0
                        powr(0) = 1d0 / r
                        do dd = 1, deg
                           pow1(dd) = xisrc1 * pow1(dd-1)
                           pow2(dd) = xisrc2 * pow2(dd-1)
                           powr(dd) = powr(dd-1) * r2inv
                        enddo

                        do term = 1, nc
                           sum = 0d0
                           do subtrm = indlo(term), indhi(term)
                              sum = sum +
     &                             cfac(subtrm) *
     &                             pow1(cpow(1, subtrm)) *
     &                             pow2(cpow(2, subtrm)) *
     &                             powr(cpow(3, subtrm))
                           enddo

!                           srcdst(term, i, j, k, patch1, patch2) =
!     &                          sum
!     petermc, 14 oct 2004:  changed sum to -sum
                           sdcrse(term, pt, patch1, patch2) = -sum

                        enddo

!     Recall that we set sum = ev(dest1, dest2) before adding new terms.
!     ev(dest1, dest2) = sum
                     enddo
                  enddo
!     end of loop over patches
               endif
!     end of check whether in range d2lo:d2hi
            enddo
         enddo
      enddo
!     end of loop over destination points

      if (pt .ne. ncpts) then
         print *, 'getcoarsematrix: ran through ', pt,
     &        ' coarse points instead of ', ncpts
         ier = 2
         return
      endif

      ier = 0

      return
      end


c=======================================================================
      subroutine getrectdirectcoarsematrix(sdcrse, wtsint,
     &     ncpts, ndstc, cadd, nsrc, hsrc, hdst, cpbase, evbase,
     &     vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     of D1 on the coarse points in a region outside of D2.
!
!     => inputs, <= outputs
!     sdcrse <= the matrix
!     nsrc => vector, lengths in cells of domain D1
!     deg => degree of multipole expansions
!     h => vector, mesh spacing in three dimensions
!     cadd => length of extra layer of coarse points beyond D2
!     cpind, cfac, cpow => from subroutine deriv1overr, describing
!                         expressions for partial derivatives of 1/r
!     nterms => number of terms in cpind, cfac, cpow
!     ncpts => number of coarse points
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer nsrc(3), ndstc(3), ncpts
      integer cadd
      real*8 wtsint(0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 sdcrse(ncpts, 0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 cpbase(3), evbase(3), hsrc(3), hdst(3)
      integer vrbos
      integer ier
!     --- locals ---
      integer idir
      integer idst, jdst, kdst
      integer isrc, jsrc, ksrc
      integer d2lo(3), d2hi(3)
      real*8 r, fac
      real*8 xsrc, ysrc, zsrc
      real*8 xdst, ydst, zdst
      real*8 xdif, ydif, zdif
      integer pt

      do idir = 1, 3
         d2lo(idir) = 0
         d2hi(idir) = ndstc(idir)
      enddo

      fac = fourth / PI

!     Loop over all coarse outer points.
      pt = 0
      do kdst = -cadd, ndstc(3) + cadd
         zdst = evbase(3) + kdst * hdst(3)
         do jdst = -cadd, ndstc(2) + cadd
            ydst = evbase(2) + jdst * hdst(2)
            do idst = -cadd, ndstc(1) + cadd
               xdst = evbase(1) + idst * hdst(1)

               if (.not. (
     &              (idst .ge. d2lo(1) .and. idst .le. d2hi(1)) .and.
     &              (jdst .ge. d2lo(2) .and. jdst .le. d2hi(2)) .and.
     &              (kdst .ge. d2lo(3) .and. kdst .le. d2hi(3)))) then

                  pt = pt + 1

!     Loop over all charges on source face.
                  do ksrc = 0, nsrc(3)
                     zsrc = cpbase(3) + ksrc * hsrc(3)
                     zdif = zdst - zsrc
                     do jsrc = 0, nsrc(2)
                        ysrc = cpbase(2) + jsrc * hsrc(2)
                        ydif = ydst - ysrc
                        do isrc = 0, nsrc(1)
                           xsrc = cpbase(1) + isrc * hsrc(1)
                           xdif = xdst - xsrc

                           r = sqrt(xdif*xdif + ydif*ydif + zdif*zdif)

                           sdcrse(pt, isrc, jsrc, ksrc) =
     &                          - wtsint(isrc, jsrc, ksrc) * fac / r
                        enddo
                     enddo
                  enddo

               endif

            enddo
         enddo
      enddo

!     End of loop over all coarse outer points.
      if (pt .ne. ncpts) then
         print *, 'getrectdirectcoarsematrix: ',
     &        'ran through ', pt,
     &        ' coarse points instead of ', ncpts
         ier = 2
         return
      endif

      ier = 0

      return
      end
