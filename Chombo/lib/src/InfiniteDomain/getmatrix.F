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

      subroutine getrectmatrix(srcdst, iev, icp, nsrc, ndst, refc,
     &     ptchsz, npchs1, npchs2,
     &     ndst1c, ndst2c, deg, layer, h,
     &     cpind, cfac, cpow, nterms, vrbos, ier)
!     This function calculates the coefficients of a matrix that
!     converts multipole expansions from patches on a face of D1
!     to their evaluation on the coarse points of a face of D2.
!
!     => input, <= outputs
!     srcdst <= the matrix
!     iev => +/- 1,2,3; destination face on which to evaluate
!     icp => +/- 1,2,3; source face with multipole expansions
!     nsrc => vector, lengths in cells of domain D1
!     ndst => vector, lengths in cells of domain D2
!     refc => coarsening ratio on the destination face
!     ptchsz => length in cells of patches
!     npchs1 => number of patches in 1st direction on face of D1
!     npchs2 => number of patches in 2nd direction on face of D1
!     ndst1c => number of coarse points in 1st direction on face of D2
!     ndst2c => number of coarse points in 2nd direction on face of D2
!     deg => degree of multipole expansions
!     layer => length of extra layer of coarse points to evaluate
!     h => vector, mesh spacing in three dimensions
!     cpind, cfac, cpow => from subroutine deriv1overr, describing
!                         expressions for partial derivatives of 1/r
!     nterms => number of terms in cpind, cfac, cpow
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer iev, icp, nsrc(3), ndst(3), refc
      integer ptchsz, npchs1, npchs2
      integer ndst1c, ndst2c, deg, layer
!     ev is at coarse points of a face of D2, plus a layer
!      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
!     cp is at centers of patches of a face of D1
!      real*8 cp(0:npchs-1, 0:npchs-1, ((deg+1)*(deg+2))/2)
      real*8 h(3)
      integer nterms
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      real*8 srcdst(((deg+1)*(deg+2))/2,
     &     -layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer idir, buf(3)
      integer patch1, patch2, dest1, dest2, nc
      integer plsrc, pldst, sidsrc, siddst
      integer isrc1, isrc2, idst1, idst2, term, subtrm
      integer endp1, endp2
      real*8 hdst1c, hdst2c
      real*8 bdst(3), dst(3)
      real*8 xisrc1, xisrc2, xplsrc, r, r2inv, sum
      integer iergpi, dd
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)
      real*8 srcpl, pctr1(0:npchs1-1), pctr2(0:npchs2-1)
      integer offlo1, offlo2

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
         print *, 'getrectmatrix:  illegal icp = ', icp
         ier = 1
         return
      endif

!     set pldst=abs(iev), siddst=sign(iev), and let idst1, idst2 be
!     the two perpendicular directions on the destination face
      call getplaneindices(iev, pldst, siddst, idst1, idst2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getrectmatrix:  illegal iev = ', iev
         ier = 2
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc, hilo(siddst), pldst
 110     format (1x, 'source face ', a2, i2, '; dest face ', a2, i2)
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

!     mesh spacing for coarsened level
      hdst1c = refc * h(idst1)
      hdst2c = refc * h(idst2)

!     physical width of each patch
!      hpch1 = ptchsz * h(isrc1)
!      hpch2 = ptchsz * h(isrc2)

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

!     We should have ndst1c == ceil(ndst(idst1) / refc).
      offlo1 = (refc * ndst1c - ndst(idst1)) / 2
      offlo2 = (refc * ndst2c - ndst(idst2)) / 2
!     Set bdst(1:3) to be physical coordinates of
!     the point with evaluation at ev(0, 0).
!     This point has coordinates (-offlo1, -offlo2)
!     in the idst1 and idst2 directions.
!     Recall that we also evaluate at points with negative
!     indices into ev:  these are points outside D2.
!     A face of D2 extends in space
!     from 0 to ndst*h in the idst1 and idst2 directions,
!     and fixed at either 0 or ndst*h in pldst direction.
      bdst(idst1) = -offlo1 * h(idst1)
      bdst(idst2) = -offlo2 * h(idst2)
      if (siddst .eq. -1) then
         bdst(pldst) = zero
      elseif (siddst .eq. 1) then
         bdst(pldst) = ndst(pldst) * h(pldst)
      endif
      if (vrbos .gt. 1) then
         print *, 'dest is based at ', bdst(1), bdst(2), bdst(3),
     &        ', fixed dimension ', pldst
      endif

!     dst(1:3) holds physical coordinates of destination point.
!     fixed dimension pldst; later we'll set
!     dst(idst1) = bdst(idst1) + dest1*hbig
!     dst(idst2) = bdst(idst2) + dest2*hbig
      dst(pldst) = bdst(pldst)

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         do patch1 = 0, npchs1-1
!     Patch with index (patch1, patch2) is centered at point with
!     coordinate pctr1(patch1) in first parallel direction,
!     coordinate pctr2(patch2) in second parallel direction,
!     coordinate srcpl in normal direction.
!     The patch has multipole expansion cp(patch1, patch2, 1:nc).

!     Loop over all destination points on destination face.
            do dest2 = -layer, ndst2c + layer
               dst(idst2) = bdst(idst2) + dest2 * hdst2c
               do dest1 = -layer, ndst1c + layer
                  dst(idst1) = bdst(idst1) + dest1 * hdst1c
!     Destination point is dst(1:3), and we add to ev(dest1, dest2).
!     xisrc1, xisrc2, xplsrc are the components of the displacement vector
!     from patch center to dst.
!     The multipole expansions are based on powers of xisrc1 and xisrc2.
                  xisrc1 = dst(isrc1) - pctr1(patch1)
                  xisrc2 = dst(isrc2) - pctr2(patch2)
                  xplsrc = dst(plsrc) - srcpl
!     r is length of x.
                  r = dsqrt(xisrc1**2 + xisrc2**2 + xplsrc**2)
                  r2inv = one / (r*r)

!     for k = 0:deg
!     pow1(k) == xisrc1**k
!     pow2(k) == xisrc2**k
!     powr(k) == 1 / r**(1+2*k)
                  pow1(0) = one
                  pow2(0) = one
                  powr(0) = one / r
                  do dd = 1, deg
                     pow1(dd) = xisrc1 * pow1(dd-1)
                     pow2(dd) = xisrc2 * pow2(dd-1)
                     powr(dd) = powr(dd-1) * r2inv
                  enddo

                  do term = 1, nc
                     sum = zero
                     do subtrm = indlo(term), indhi(term)
                        sum = sum +
     &                    cfac(subtrm) *
     &                    pow1(cpow(1, subtrm)) *
     &                    pow2(cpow(2, subtrm)) *
     &                    powr(cpow(3, subtrm))
                     enddo
                     srcdst(term, dest1, dest2, patch1, patch2) = sum
                  enddo

!                     newtrm =
!     &                    cp(patch1, patch2, cpind(term)) *
!     &                    cfac(term) *
!     &                    (xisrc1**cpow(1, term)) *
!     &                    (xisrc2**cpow(2, term)) /
!     &                    (r**(1 + 2*cpow(3, term)))

!     Recall that we set sum = ev(dest1, dest2) before adding new terms.
!                  ev(dest1, dest2) = sum

               enddo
            enddo
!     end of loop over destination points
         enddo
      enddo
!     end of loop over patches

      ier = 0

      return
      end


!=======================================================================
      subroutine getplaneindices(ind, plane, side, i1, i2, ier)
!     => inputs, <= outputs
!     ind => +/- 1, 2, 3
!     plane <= which dimension is constant (1, 2, 3)
!     side <= -1 for low, +1 for high
!     i1, i2 <= the two variable dimensions (1, 2, 3 but not plane)
!     ier <= return code, 0 if and only if success
      implicit none
      integer ind, plane, side, i1, i2, ier
      integer iabs

      plane = iabs(ind)
      side = ind / plane

      if (side .ne. 1 .and. side .ne. -1) then
!     illegal side
         ier = 1
         return
      endif

      if (plane .eq. 1) then
!     plane x = const
         i1 = 2
         i2 = 3
      elseif (plane .eq. 2) then
!     plane y = const
         i1 = 1
         i2 = 3
      elseif (plane .eq. 3) then
!     plane z = const
         i1 = 1
         i2 = 2
      else
!     illegal plane
         ier = 2
         return
      endif

      ier = 0

      return
      end


!=======================================================================
      subroutine getcenters(pctr, npchs, bsrc, ptchsz, endpsz, h)
!     This subroutine returns the physical coordinates of the centers of
!     patches in one dimension.
!
!     => inputs, <= outputs
!     pctr(0:npchs-1) <= centers of patches in this dimension
!     npchs => number of patches in this dimension
!     bsrc => base physical coordinate of patch in this dimension
!     ptchsz => regular patch length, for patches 1:npchs-2
!     endpsz => length of patches 0 and npchs-1
!     h => mesh spacing in this dimension
      implicit none
!     --- arguments ---
      integer npchs, ptchsz, endpsz
      real*8 pctr(0:npchs-1), bsrc, h
!     --- locals ---
      integer patch
      real*8 hpatch, pchoff

      pchoff = endpsz * h
!     first (index 0) patch extends
!     from bsrc
!     to bsrc + pchoff
      pctr(0) = bsrc + pchoff * half
      if (npchs .gt. 1) then
!     hpatch = length of regular patch
         hpatch = ptchsz * h
         do patch = 1, npchs-2
!     regular patch extends
!     from bsrc + pchoff + (patch-1)*hpatch
!     to bsrc + pchoff + patch*hpatch
            pctr(patch) = bsrc + pchoff + patch * hpatch - half * hpatch
         enddo
!     last (index npchs-1) patch extends
!     from bsrc + pchoff + (npchs-2)*hpatch
!     to bsrc + pchoff + (npchs-2)*hpatch + pchoff
         pctr(npchs-1) = bsrc + pchoff*three*half + (npchs-2)*hpatch
      endif

      return
      end


!=======================================================================
      subroutine getdirectrectmatrix(srcdst, wtsint, iev, icp, nsrc,
     &     ndst1c, ndst2c, layer, hsrc, hdst, cpbase, evbase,
     &     vrbos, ier)
!     This function calculates the coefficients of a matrix that
!     converts multipole expansions from patches on a face of D1
!     to their evaluation on the coarse points of a face of D2.
!
!     => input, <= outputs
!     srcdst <= the matrix
!     wtsint => integration weights on charges
!     iev => +/- 1,2,3; destination face on which to evaluate
!     icp => +/- 1,2,3; source face with multipole expansions
!     nsrc => vector, lengths in cells of domain D1
!     ndst1c => number of coarse points in 1st direction on face of D2
!     ndst2c => number of coarse points in 2nd direction on face of D2
!     layer => length of extra layer of coarse points to evaluate
!     hsrc => mesh spacings on source face D1
!     hdst => mesh spacings on destination domain D2
!     cpbase => physical coordinate of low corner of source face
!     evbase => physical coordinate of low corner of destination face
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer iev, icp, nsrc(3)
      integer ndst1c, ndst2c, layer
!     ev is at coarse points of a face of D2, plus a layer
!      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
      real*8 wtsint(0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 srcdst(-layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 cpbase(3), evbase(3), hsrc(3), hdst(3)
      integer vrbos
      integer ier
!     --- locals ---
      integer src1, src2, src3, dest1, dest2
      integer plsrc, pldst, sidsrc, siddst
      integer isrc1, isrc2, idst1, idst2
      real*8 dst(3), src(3)
      real*8 fac, dif1, dif2, dif3, r, wthere
      integer iergpi
      character*2 hilo(-1:1)

!     set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
!     the two perpendicular directions on the source face
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getrectmatrix:  illegal icp = ', icp
         ier = 1
         return
      endif

!     set pldst=abs(iev), siddst=sign(iev), and let idst1, idst2 be
!     the two perpendicular directions on the destination face
      call getplaneindices(iev, pldst, siddst, idst1, idst2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getrectmatrix:  illegal iev = ', iev
         ier = 2
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc, hilo(siddst), pldst
 110     format (1x, 'source face ', a2, i2, '; dest face ', a2, i2)
      endif

!     dst(1:3) holds physical coordinates of destination point.
!     fixed dimension pldst; later we'll set
!     dst(idst1) = evbase(idst1) + dest1 * hdst(idst1)
!     dst(idst2) = evbase(idst2) + dest2 * hdst(idst2)
      dst(pldst) = evbase(pldst)

      fac = fourth / PI

!     Loop over all charges on source face.
      do src3 = 0, nsrc(3)
         src(3) = cpbase(3) + src3 * hsrc(3)
         do src2 = 0, nsrc(2)
            src(2) = cpbase(2) + src2 * hsrc(2)
            do src1 = 0, nsrc(1)
               src(1) = cpbase(1) + src1 * hsrc(1)

               wthere = wtsint(src1, src2, src3)

!     Loop over all coarse points on destination face.
               do dest2 = -layer, ndst2c + layer
                  dst(idst2) = evbase(idst2) + dest2 * hdst(idst2)
                  do dest1 = -layer, ndst1c + layer
                     dst(idst1) = evbase(idst1) + dest1 * hdst(idst1)

                     dif1 = dst(1) - src(1)
                     dif2 = dst(2) - src(2)
                     dif3 = dst(3) - src(3)

                     r = sqrt(dif1**2 + dif2**2 + dif3**2)

                     srcdst(dest1, dest2, src1, src2, src3) =
     &                    wthere * fac / r

                  enddo
               enddo
!     end of loop over destination points
            enddo
         enddo
      enddo
!     end of loop over source points

      ier = 0

      return
      end

!=======================================================================
!     For debugging.
      subroutine getorirectmatrix(ori,
     &     srcdst, iev, icp, nsrc, ndst, refc,
     &     ptchsz, npchs1, npchs2,
     &     ndst1c, ndst2c, deg, layer, h,
     &     cpind, cfac, cpow, nterms, vrbos, ier)
!     This function calculates the coefficients of a matrix that
!     converts multipole expansions from patches on a face of D1
!     to their evaluation on the coarse points of a face of D2.
!
!     => input, <= outputs
!     srcdst <= the matrix
!     iev => +/- 1,2,3; destination face on which to evaluate
!     icp => +/- 1,2,3; source face with multipole expansions
!     nsrc => vector, lengths in cells of domain D1
!     ndst => vector, lengths in cells of domain D2
!     refc => coarsening ratio on the destination face
!     ptchsz => length in cells of patches
!     npchs1 => number of patches in 1st direction on face of D1
!     npchs2 => number of patches in 2nd direction on face of D1
!     ndst1c => number of coarse points in 1st direction on face of D2
!     ndst2c => number of coarse points in 2nd direction on face of D2
!     deg => degree of multipole expansions
!     layer => length of extra layer of coarse points to evaluate
!     h => vector, mesh spacing in three dimensions
!     cpind, cfac, cpow => from subroutine deriv1overr, describing
!                         expressions for partial derivatives of 1/r
!     nterms => number of terms in cpind, cfac, cpow
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer ori, iev, icp, nsrc(3), ndst(3), refc
      integer ptchsz, npchs1, npchs2
      integer ndst1c, ndst2c, deg, layer
!     ev is at coarse points of a face of D2, plus a layer
!      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
!     cp is at centers of patches of a face of D1
!      real*8 cp(0:npchs-1, 0:npchs-1, ((deg+1)*(deg+2))/2)
      real*8 h(3)
      integer nterms
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      real*8 srcdst(((deg+1)*(deg+2))/2,
     &     -layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer idir, buf(3)
      integer patch1, patch2, dest1, dest2, nc
      integer plsrc, pldst, sidsrc, siddst
      integer isrc1, isrc2, idst1, idst2, term, subtrm
      integer endp1, endp2
      real*8 hdst1c, hdst2c
      real*8 bdst(3), dst(3)
      real*8 xisrc1, xisrc2, xplsrc, r, r2inv, sum
      integer iergpi, dd
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)
      real*8 srcpl, pctr1(0:npchs1-1), pctr2(0:npchs2-1)
      integer offlo1, offlo2, t1, t2, sgn1, sgn2
      integer flip1(((deg+1)*(deg+2))/2)
      integer flip2(((deg+1)*(deg+2))/2)
      integer flip12(((deg+1)*(deg+2))/2)
      integer flipme(((deg+1)*(deg+2))/2)
      integer p1beg, p1end, p1inc
      integer p2beg, p2end, p2inc
      integer d1beg, d1end, d1inc
      integer d2beg, d2end, d2inc

!     original, dest flip, source flip.
!     flip is 0 if none, 1 if first coordinate, 2 if second coordinate.
      integer ors(0:2, 0:35)
      data ors /
     &     0, 0, 0,
     &     1, 0, 0,
     &     2, 0, 0,
     &     18, 0, 0,
     &     1, 1, 0,
     &     2, 2, 0,
     &     6, 0, 0,
     &     7, 0, 0,
     &     8, 0, 0,
     &     6, 1, 0,
     &     25, 0, 0,
     &     8, 2, 0,
     &     12, 0, 0,
     &     13, 0, 0,
     &     14, 0, 0,
     &     12, 1, 0,
     &     13, 2, 0,
     &     32, 0, 0,
     &     18, 0, 0,
     &     1, 0, 1,
     &     2, 0, 1,
     &     0, 0, 0,
     &     1, 1, 1,
     &     2, 2, 1,
     &     6, 0, 1,
     &     25, 0, 0,
     &     8, 0, 2,
     &     6, 1, 1,
     &     7, 0, 0,
     &     8, 2, 2,
     &     12, 0, 2,
     &     13, 0, 2,
     &     32, 0, 0,
     &     12, 1, 2,
     &     13, 2, 2,
     &     14, 0, 0 /

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
         print *, 'getrectmatrix:  illegal icp = ', icp
         ier = 1
         return
      endif

!     set pldst=abs(iev), siddst=sign(iev), and let idst1, idst2 be
!     the two perpendicular directions on the destination face
      call getplaneindices(iev, pldst, siddst, idst1, idst2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getrectmatrix:  illegal iev = ', iev
         ier = 2
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc, hilo(siddst), pldst
 110     format (1x, 'source face ', a2, i2, '; dest face ', a2, i2)
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

!     mesh spacing for coarsened level
      hdst1c = refc * h(idst1)
      hdst2c = refc * h(idst2)

!     physical width of each patch
!      hpch1 = ptchsz * h(isrc1)
!      hpch2 = ptchsz * h(isrc2)

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

!     We should have ndst1c == ceil(ndst(idst1) / refc).
      offlo1 = (refc * ndst1c - ndst(idst1)) / 2
      offlo2 = (refc * ndst2c - ndst(idst2)) / 2
!     Set bdst(1:3) to be physical coordinates of
!     the point with evaluation at ev(0, 0).
!     This point has coordinates (-offlo1, -offlo2)
!     in the idst1 and idst2 directions.
!     Recall that we also evaluate at points with negative
!     indices into ev:  these are points outside D2.
!     A face of D2 extends in space
!     from 0 to ndst*h in the idst1 and idst2 directions,
!     and fixed at either 0 or ndst*h in pldst direction.
      bdst(idst1) = -offlo1 * h(idst1)
      bdst(idst2) = -offlo2 * h(idst2)
      if (siddst .eq. -1) then
         bdst(pldst) = zero
      elseif (siddst .eq. 1) then
         bdst(pldst) = ndst(pldst) * h(pldst)
      endif
      if (vrbos .gt. 1) then
         print *, 'dest is based at ', bdst(1), bdst(2), bdst(3),
     &        ', fixed dimension ', pldst
      endif

!     dst(1:3) holds physical coordinates of destination point.
!     fixed dimension pldst; later we'll set
!     dst(idst1) = bdst(idst1) + dest1*hbig
!     dst(idst2) = bdst(idst2) + dest2*hbig
      dst(pldst) = bdst(pldst)

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         do patch1 = 0, npchs1-1
!     Patch with index (patch1, patch2) is centered at point with
!     coordinate pctr1(patch1) in first parallel direction,
!     coordinate pctr2(patch2) in second parallel direction,
!     coordinate srcpl in normal direction.
!     The patch has multipole expansion cp(patch1, patch2, 1:nc).

!     Loop over all destination points on destination face.
            do dest2 = -layer, ndst2c + layer
               dst(idst2) = bdst(idst2) + dest2 * hdst2c
               do dest1 = -layer, ndst1c + layer
                  dst(idst1) = bdst(idst1) + dest1 * hdst1c
!     Destination point is dst(1:3), and we add to ev(dest1, dest2).
!     xisrc1, xisrc2, xplsrc are the components of the displacement vector
!     from patch center to dst.
!     The multipole expansions are based on powers of xisrc1 and xisrc2.
                  xisrc1 = dst(isrc1) - pctr1(patch1)
                  xisrc2 = dst(isrc2) - pctr2(patch2)
                  xplsrc = dst(plsrc) - srcpl
!     r is length of x.
                  r = dsqrt(xisrc1**2 + xisrc2**2 + xplsrc**2)
                  r2inv = one / (r*r)

!     for k = 0:deg
!     pow1(k) == xisrc1**k
!     pow2(k) == xisrc2**k
!     powr(k) == 1 / r**(1+2*k)
                  pow1(0) = one
                  pow2(0) = one
                  powr(0) = one / r
                  do dd = 1, deg
                     pow1(dd) = xisrc1 * pow1(dd-1)
                     pow2(dd) = xisrc2 * pow2(dd-1)
                     powr(dd) = powr(dd-1) * r2inv
                  enddo

                  do term = 1, nc
                     sum = zero
                     do subtrm = indlo(term), indhi(term)
                        sum = sum +
     &                    cfac(subtrm) *
     &                    pow1(cpow(1, subtrm)) *
     &                    pow2(cpow(2, subtrm)) *
     &                    powr(cpow(3, subtrm))
                     enddo
                     srcdst(term, dest1, dest2, patch1, patch2) = sum
                  enddo

!                     newtrm =
!     &                    cp(patch1, patch2, cpind(term)) *
!     &                    cfac(term) *
!     &                    (xisrc1**cpow(1, term)) *
!     &                    (xisrc2**cpow(2, term)) /
!     &                    (r**(1 + 2*cpow(3, term)))

!     Recall that we set sum = ev(dest1, dest2) before adding new terms.
!                  ev(dest1, dest2) = sum

               enddo
            enddo
!     end of loop over destination points
         enddo
      enddo
!     end of loop over patches

      term = 0
      do t1 = 0, deg
         if (mod(t1, 2) .eq. 0) then
            sgn1 = 1
         else
            sgn1 = -1
         endif
         do t2 = 0, deg
            if (t1 + t2 .le. deg) then
               term = term + 1
               if (mod(t2, 2) .eq. 0) then
                  sgn2 = 1
               else
                  sgn2 = -1
               endif
               flip1(term) = sgn1
               flip2(term) = sgn2
               flip12(term) = sgn1 * sgn2
            endif
         enddo
      enddo


      p1beg = 0
      p1end = npchs1-1
      p1inc = 1
      p2beg = 0
      p2end = npchs2-1
      p2inc = 1
      if (ors(2, ori) .eq. 0) then
         do term = 1, nc
            flipme(term) = 1
         enddo
      elseif (ors(2, ori) .eq. 1) then
         p1beg = npchs1-1
         p1end = 0
         p1inc = -1
         do term = 1, nc
            flipme(term) = flip1(term)
         enddo
      elseif (ors(2, ori) .eq. 2) then
         p2beg = npchs2-1
         p2end = 0
         p2inc = -1
         do term = 1, nc
            flipme(term) = flip2(term)
         enddo
      else
         print *, 'error in ors'
         ier = -1
         return
      endif

      d1beg = -layer
      d1inc = 1
      d1end = ndst1c + layer
      d2beg = -layer
      d2end = ndst2c + layer
      d2inc = 1
      if (ors(1, ori) .eq. 0) then
c     nothing
      elseif (ors(1, ori) .eq. 1) then
         d1beg = ndst1c + layer
         d1end = -layer
         d1inc = -1
      elseif (ors(1, ori) .eq. 2) then
         d2beg = ndst2c + layer
         d2end = -layer
         d2inc = -1
      else
         print *, 'error in ors'
         ier = -1
         return
      endif

      if (.false.) then
      do patch2 = p2beg, p2end, p2inc
         do patch1 = p1beg, p1end, p1inc
            do dest2 = d2beg, d2end, d2inc
               do dest1 = d1beg, d1end, d1inc
                  do term = 1, nc
                     print *, 'ori ', ori, flipme(term) *
     &                    srcdst(term, dest1, dest2, patch1, patch2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      endif

      ier = 0

      return
      end
