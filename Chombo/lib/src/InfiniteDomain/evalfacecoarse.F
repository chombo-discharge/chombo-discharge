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

      subroutine evalfacecoarse(ev, cp,
     &     npchs, n2c, nc, layer, ori, srcdst, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     on the coarse points of another face, and adds the evaluation to ev.
!
!     => inputs, <= outputs
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     npchs => number of patches in each dimension
!     n2c => ceil(n2/ptchsz), number of coarse points in each dimension
!     nc => number of coefficients in multipole expansions
!     layer => length of extra layer of coarse points to evaluate
!     ori => the particular orientation of source and destination
!     srcdst => input describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer npchs, n2c, nc, layer, ori
!     ev is at coarse points of a face of D2, plus a layer
      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
!     cp is at centers of patches of a face of D1
!     petermc, 5 Jan 2004:  new order (coefficient, patch1, patch2)
!     petermc, 26 Aug 2004:  new order (patch1, patch2, coefficient)
      real*8 cp(0:npchs-1, 0:npchs-1, nc)
      real*8 srcdst(nc,
     &     -layer:n2c+layer, -layer:n2c+layer,
     &     0:npchs-1, 0:npchs-1, 36)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, dest1, dest2, term
      real*8 newtrm

!     number of coefficients
!      nc = ((deg+1)*(deg+2))/2
!     length of D
!      n = n1 - 2*s1
!     length of D2
!      n2 = n1 + 2*s2
!     mesh spacing for coarsened level
!      hbig = ptchsz * h

!     Loop over all patches on source face.
      do patch2 = 0, npchs-1
         do patch1 = 0, npchs-1
!     Loop over all coarse points on destination face.
            do dest2 = -layer, n2c + layer
               do dest1 = -layer, n2c + layer
!     Loop over all terms.
                  do term = 1, nc
                     newtrm = cp(patch1, patch2, term) *
     &                    srcdst(term,
     &                    dest1, dest2,
     &                    patch1, patch2, ori)
                     ev(dest1, dest2) = ev(dest1, dest2) + newtrm
                  enddo
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

      subroutine addeval1facecoarse(ev, cp,
     &     npchs1, npchs2, ndst1c, ndst2c,
     &     nc, layer, srcdst, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     on the coarse points of another face, and adds the evaluation to ev.
!
!     => inputs, <= outputs
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     npchs1 => number of patches in first source-face dimension
!     npchs2 => number of patches in second source-face dimension
!     ndst1c => number of coarse points in first dest-face dimension
!     ndst2c => number of coarse points in second dest-face dimension
!     nc => number of coefficients in multipole expansions
!     layer => length of extra layer of coarse points to evaluate
!     srcdst => input describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer npchs1, npchs2, ndst1c, ndst2c, nc, layer
!     ev is at coarse points of a face of D2, plus a layer
      real*8 ev(-layer:ndst1c+layer, -layer:ndst2c+layer)
!     cp is at centers of patches of a face of D1
!     petermc, 5 Jan 2004:  new order (coefficient, patch1, patch2)
!     petermc, 26 Aug 2004:  new order (patch1, patch2, coefficient)
      real*8 cp(0:npchs1-1, 0:npchs2-1, nc)
      real*8 srcdst(nc,
     &     -layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, dest1, dest2, term
      real*8 newtrm

!     number of coefficients
!      nc = ((deg+1)*(deg+2))/2
!     length of D
!      n = n1 - 2*s1
!     length of D2
!      n2 = n1 + 2*s2
!     mesh spacing for coarsened level
!      hbig = ptchsz * h

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         do patch1 = 0, npchs1-1
!     Loop over all coarse points on destination face.
            do dest2 = -layer, ndst2c + layer
               do dest1 = -layer, ndst1c + layer
!     Loop over all terms.
                  do term = 1, nc
                     newtrm = cp(patch1, patch2, term) *
     &                    srcdst(term, dest1, dest2, patch1, patch2)
                     ev(dest1, dest2) = ev(dest1, dest2) + newtrm
                  enddo
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

      subroutine addevalsym1facecoarse(ev, cp, dflip, sflip,
     &     npchs1, npchs2, ndst1c, ndst2c,
     &     deg, nc, layer, srcdst, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     on the coarse points of another face, and adds the evaluation to ev.
!
!     => inputs, <= outputs
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     dflip => 1 if flip first dest coordinate, 2 if second, 0 if neither
!     sflip => 1 if flip first source coordinate, 2 if second, 0 if neither
!     npchs1 => number of patches in first source-face dimension
!     npchs2 => number of patches in second source-face dimension
!     ndst1c => number of coarse points in first dest-face dimension
!     ndst2c => number of coarse points in second dest-face dimension
!     deg => degree of multipole expansions
!     nc => number of coefficients in multipole expansions
!     layer => length of extra layer of coarse points to evaluate
!     srcdst => input describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer dflip, sflip
      integer npchs1, npchs2, ndst1c, ndst2c, deg, nc, layer
!     ev is at coarse points of a face of D2, plus a layer
      real*8 ev(-layer:ndst1c+layer, -layer:ndst2c+layer)
!     cp is at centers of patches of a face of D1
!     petermc, 5 Jan 2004:  new order (coefficient, patch1, patch2)
!     petermc, 26 Aug 2004:  new order (patch1, patch2, coefficient)
      real*8 cp(0:npchs1-1, 0:npchs2-1, nc)
      real*8 srcdst(nc,
     &     -layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, dest1, dest2, term
      integer p1mat, p2mat, d1mat, d2mat
      real*8 newtrm
      integer t1, t2, sgntrm, sgnmod(0:1), flipme(nc)
      data sgnmod / 1, -1 /

!     number of coefficients
!      nc = ((deg+1)*(deg+2))/2
!     length of D
!      n = n1 - 2*s1
!     length of D2
!      n2 = n1 + 2*s2
!     mesh spacing for coarsened level
!      hbig = ptchsz * h

      if (sflip .eq. 0) then
         do term = 1, nc
            flipme(term) = 1
         enddo
      elseif (sflip .eq. 1) then
         term = 0
         do t1 = 0, deg
!     If t1 is even then mod == 0 and sgntrm = 1.
!     If t1 is odd then mod == 1 and sgntrm = -1.
            sgntrm = sgnmod(mod(t1, 2))
            do t2 = 0, deg
               if (t1 + t2 .le. deg) then
                  term = term + 1
                  flipme(term) = sgntrm
               endif
            enddo
         enddo
      elseif (sflip .eq. 2) then
         term = 0
         do t1 = 0, deg
            do t2 = 0, deg
               if (t1 + t2 .le. deg) then
!     If t2 is even then mod == 0 and sgntrm = 1.
!     If t2 is odd then mod == 1 and sgntrm = -1.
                  sgntrm = sgnmod(mod(t2, 2))
                  term = term + 1
                  flipme(term) = sgntrm
               endif
            enddo
         enddo
      else
         print *, 'addevalsym1facecoarse:  sflip must be 0 or 1 or 2'
         ier = -1
         return
      endif

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         if (sflip .eq. 2) then
            p2mat = npchs2-1 - patch2
         else
            p2mat = patch2
         endif
         do patch1 = 0, npchs1-1
            if (sflip .eq. 1) then
               p1mat = npchs1-1 - patch1
            else
               p1mat = patch1
            endif
!     Loop over all coarse points on destination face.
            do dest2 = -layer, ndst2c + layer
               if (dflip .eq. 2) then
                  d2mat = ndst2c - dest2
               else
                  d2mat = dest2
               endif
               do dest1 = -layer, ndst1c + layer
                  if (dflip .eq. 1) then
                     d1mat = ndst1c - dest1
                  else
                     d1mat = dest1
                  endif
!     Loop over all terms.
                  do term = 1, nc
                     newtrm = cp(patch1, patch2, term) *
     &                    flipme(term) *
     &                    srcdst(term, d1mat, d2mat, p1mat, p2mat)
                     ev(dest1, dest2) = ev(dest1, dest2) + newtrm
                  enddo
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

      subroutine addeval1directfacecoarse(ev, cp,
     &     nsrc, ndst1c, ndst2c, layer, srcdst,
     &     vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     on the coarse points of another face, and adds the evaluation to ev.
!
!     => inputs, <= outputs
!     ev <= evaluation
!     cp => charges on a face
!     ndst1c => number of coarse points in first dest-face dimension
!     ndst2c => number of coarse points in second dest-face dimension
!     layer => length of extra layer of coarse points to evaluate
!     srcdst => input describing coefficients 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer nsrc(3), ndst1c, ndst2c, layer
!     ev is at coarse points of a face of D2, plus a layer
      real*8 ev(-layer:ndst1c+layer, -layer:ndst2c+layer)
      real*8 cp(0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 srcdst(-layer:ndst1c+layer, -layer:ndst2c+layer,
     &     0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      integer vrbos
      integer ier
!     --- locals ---
      integer src1, src2, src3, dest1, dest2

!     Loop over all charges on source face.
      do src3 = 0, nsrc(3)
         do src2 = 0, nsrc(2)
            do src1 = 0, nsrc(1)
!     Loop over all coarse points on destination face.
               do dest2 = -layer, ndst2c + layer
                  do dest1 = -layer, ndst1c + layer
                     ev(dest1, dest2) = ev(dest1, dest2) +
     &                    cp(src1, src2, src3) *
     &                    srcdst(dest1, dest2, src1, src2, src3)
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
