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

      subroutine evaloutercoarse(ev, cp,
     &     ncpts, npchs, n2c, cadd, nc, ori, sdcrse, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     on the coarse points of another face, and adds the evaluation to ev.
!
!     => inputs, <= outputs
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     ncpts => number of coarse outer points
!     npchs => number of patches in each dimension
!     n2c => floor(n2/refc), coarse length of D2
!     cadd => width of extra layer of coarse points beyond D2
!     nc => number of coefficients in multipole expansions
!     ori => the particular orientation of source and destination
!     srcdst => array describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer ncpts, npchs, n2c, cadd, nc, ori
!     ev is at coarse outer points beyond D2
      real*8 ev(-cadd:n2c+cadd, -cadd:n2c+cadd, -cadd:n2c+cadd)
!     cp is at centers of patches of a face of D1
      real*8 cp(0:npchs-1, 0:npchs-1, nc)
      real*8 sdcrse(nc, ncpts, 0:npchs-1, 0:npchs-1, 6)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, d2lo, d2hi, pt, term, i, j, k
      real*8 newtrm

!     This is much slower if you nest the loops over source-face patches
!     within the loops over coarse destination points.

      d2lo = 0
      d2hi = n2c

!     First zero out the coarse outer points.
      do k = -cadd, n2c + cadd
         do j = -cadd, n2c + cadd
            do i = -cadd, n2c + cadd
               if (.not. (
     &              (i .ge. d2lo .and. i .le. d2hi) .and.
     &              (j .ge. d2lo .and. j .le. d2hi) .and.
     &              (k .ge. d2lo .and. k .le. d2hi))) then
                  ev(i, j, k) = zero
               endif
            enddo
         enddo
      enddo

!     Loop over all patches on source face.
      do patch2 = 0, npchs-1
         do patch1 = 0, npchs-1
!     Loop over all coarse outer points.
            pt = 0
            do k = -cadd, n2c + cadd
               do j = -cadd, n2c + cadd
                  do i = -cadd, n2c + cadd
                     if (.not. (
     &                    (i .ge. d2lo .and. i .le. d2hi) .and.
     &                    (j .ge. d2lo .and. j .le. d2hi) .and.
     &                    (k .ge. d2lo .and. k .le. d2hi))) then

                        pt = pt + 1
!                        print *, 'pt', pt, ' is ', i, j, k
!     Loop over all terms.
                        do term = 1, nc
                           newtrm = cp(patch1, patch2, term) *
     &                          sdcrse(term, pt, patch1, patch2, ori)
                           ev(i, j, k) = ev(i, j, k) + newtrm
                        enddo
                     endif

                  enddo
               enddo
            enddo
!     End of loop over all coarse outer points.
            if (pt .ne. ncpts) then
               print *, 'evaloutercoarse: ran through ', pt,
     &              ' coarse points instead of ', ncpts
               ier = 2
               return
            endif

         enddo
      enddo
!     End of loop over patches.

      ier = 0

      return
      end

c=======================================================================

      subroutine evalrectoutercoarse(ev, cp,
     &     ncpts, npchs1, npchs2, ndstc, cadd, nc,
     &     sdcrse, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     at outer coarse points.
!
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     ncpts => number of coarse outer points
!     npchs1 => number of patches in 1st direction on face of D1
!     npchs2 => number of patches in 2nd direction on face of D1
!     ndstc => floor(ndst/refc), coarse length of D2
!     cadd => width of extra layer of coarse points beyond D2
!     nc => number of coefficients in multipole expansions
!     sdcrse => array describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer ncpts, npchs1, npchs2, ndstc(3), cadd, nc
!     ev is at coarse outer points beyond D2
      real*8 ev(-cadd:ndstc(1)+cadd,
     &     -cadd:ndstc(2)+cadd,
     &     -cadd:ndstc(3)+cadd)
!     cp is at centers of patches of a face of D1
      real*8 cp(0:npchs1-1, 0:npchs2-1, nc)
      real*8 sdcrse(nc, ncpts, 0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, d2lo(3), d2hi(3), pt, term, i, j, k
      integer idir
      real*8 newtrm

!     This is much slower if you nest the loops over source-face patches
!     within the loops over coarse destination points.

      do idir = 1, 3
         d2lo(idir) = 0
         d2hi(idir) = ndstc(idir)
      enddo

!     First zero out the coarse outer points.
      do k = -cadd, ndstc(3) + cadd
         do j = -cadd, ndstc(2) + cadd
            do i = -cadd, ndstc(1) + cadd
               if (.not. (
     &              (i .ge. d2lo(1) .and. i .le. d2hi(1)) .and.
     &              (j .ge. d2lo(2) .and. j .le. d2hi(2)) .and.
     &              (k .ge. d2lo(3) .and. k .le. d2hi(3)))) then
                  ev(i, j, k) = zero
               endif
            enddo
         enddo
      enddo

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         do patch1 = 0, npchs1-1
!     Loop over all coarse outer points.
            pt = 0
            do k = -cadd, ndstc(3) + cadd
               do j = -cadd, ndstc(2) + cadd
                  do i = -cadd, ndstc(1) + cadd

                     if (.not. (
     &                    (i .ge. d2lo(1) .and. i .le. d2hi(1)) .and.
     &                    (j .ge. d2lo(2) .and. j .le. d2hi(2)) .and.
     &                    (k .ge. d2lo(3) .and. k .le. d2hi(3)))) then
                        pt = pt + 1

!     Loop over all terms.
                        do term = 1, nc
                           newtrm = cp(patch1, patch2, term) *
     &                          sdcrse(term, pt, patch1, patch2)
                           ev(i, j, k) = ev(i, j, k) + newtrm
                        enddo
                     endif

                  enddo
               enddo
            enddo
!     End of loop over all coarse outer points.
            if (pt .ne. ncpts) then
               print *, 'evalrectoutercoarse: ran through ', pt,
     &              ' coarse points instead of ', ncpts
               ier = 2
               return
            endif

         enddo
      enddo
!     End of loop over patches.

      ier = 0

      return
      end


c=======================================================================

      subroutine evalsymrectoutercoarse(ev, cp, dflip,
     &     ncpts, npchs1, npchs2, ndstc, cadd, nc,
     &     sdcrse, vrbos, ier)
!     This function evaluates multipole expansions from patches on a face
!     at outer coarse points.
!
!     ev <= evaluation
!     cp => multipole expansions on patches of a face
!     dflip => 1 if flip x, 2 if flip y, 3 if flip z, 0 if no flipping
!     ncpts => number of coarse outer points
!     npchs1 => number of patches in 1st direction on face of D1
!     npchs2 => number of patches in 2nd direction on face of D1
!     ndstc => floor(ndst/refc), coarse length of D2
!     cadd => width of extra layer of coarse points beyond D2
!     nc => number of coefficients in multipole expansions
!     sdcrse => array describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer dflip, ncpts, npchs1, npchs2, ndstc(3), cadd, nc
!     ev is at coarse outer points beyond D2
      real*8 ev(-cadd:ndstc(1)+cadd,
     &     -cadd:ndstc(2)+cadd,
     &     -cadd:ndstc(3)+cadd)
!     cp is at centers of patches of a face of D1
      real*8 cp(0:npchs1-1, 0:npchs2-1, nc)
      real*8 sdcrse(nc, ncpts, 0:npchs1-1, 0:npchs2-1)
      integer vrbos
      integer ier
!     --- locals ---
      integer patch1, patch2, d2lo(3), d2hi(3), pt, term, i, j, k
      integer idir
      integer ifirst, ilast, istep
      integer jfirst, jlast, jstep
      integer kfirst, klast, kstep
      real*8 newtrm

!     This is much slower if you nest the loops over source-face patches
!     within the loops over coarse destination points.

      do idir = 1, 3
         d2lo(idir) = 0
         d2hi(idir) = ndstc(idir)
      enddo

!     First zero out the coarse outer points.
      do k = -cadd, ndstc(3) + cadd
         do j = -cadd, ndstc(2) + cadd
            do i = -cadd, ndstc(1) + cadd
               if (.not. (
     &              (i .ge. d2lo(1) .and. i .le. d2hi(1)) .and.
     &              (j .ge. d2lo(2) .and. j .le. d2hi(2)) .and.
     &              (k .ge. d2lo(3) .and. k .le. d2hi(3)))) then
                  ev(i, j, k) = zero
               endif
            enddo
         enddo
      enddo

      ifirst = -cadd
      ilast = ndstc(1) + cadd
      istep = 1
      jfirst = -cadd
      jlast = ndstc(2) + cadd
      jstep = 1
      kfirst = -cadd
      klast = ndstc(3) + cadd
      kstep = 1
      if (dflip .eq. 1) then
         ifirst = ndstc(1) + cadd
         ilast = -cadd
         istep = -1
      elseif (dflip .eq. 2) then
         jfirst = ndstc(2) + cadd
         jlast = -cadd
         jstep = -1
      elseif (dflip .eq. 3) then
         kfirst = ndstc(3) + cadd
         klast = -cadd
         kstep = -1
!     dflip == 0 is only other legal possibility
      elseif (dflip .ne. 0) then
         print *, 'evalsymrectoutercoarse: illegal dflip = ', dflip
         ier = 3
      endif

!     Loop over all patches on source face.
      do patch2 = 0, npchs2-1
         do patch1 = 0, npchs1-1
!     Loop over all coarse outer points.
!     pt is index into sdcrse
            pt = 0
!            do k = -cadd, ndstc(3) + cadd
!               do j = -cadd, ndstc(2) + cadd
!                  do i = -cadd, ndstc(1) + cadd
            do k = kfirst, klast, kstep
               do j = jfirst, jlast, jstep
                  do i = ifirst, ilast, istep

                     if (.not. (
     &                    (i .ge. d2lo(1) .and. i .le. d2hi(1)) .and.
     &                    (j .ge. d2lo(2) .and. j .le. d2hi(2)) .and.
     &                    (k .ge. d2lo(3) .and. k .le. d2hi(3)))) then
                        pt = pt + 1

!     Loop over all terms.
                        do term = 1, nc
                           newtrm = cp(patch1, patch2, term) *
     &                          sdcrse(term, pt, patch1, patch2)
                           ev(i, j, k) = ev(i, j, k) + newtrm
                        enddo
                     endif

                  enddo
               enddo
            enddo
!     End of loop over all coarse outer points.
            if (pt .ne. ncpts) then
               print *, 'evalsymrectoutercoarse: ran through ', pt,
     &              ' coarse points instead of ', ncpts
               ier = 2
               return
            endif

         enddo
      enddo
!     End of loop over patches.

      ier = 0

      return
      end


c=======================================================================

      subroutine evaldirectrectoutercoarse(ev, cp,
     &     ncpts, ndstc, cadd, nsrc, sdcrse,
     &     vrbos, ier)
!     This function evaluates directly the potential at outer
!     coarse points due to charges on a face.
!
!     ev <= evaluation
!     cp => weighted charges on a face
!     ncpts => number of coarse outer points
!     ndstc => floor(ndst/refc), coarse length of D2
!     cadd => width of extra layer of coarse points beyond D2
!     nsrc => length in cells of face of D1
!     sdcrse => array describing expansions for partial derivatives of 1/r
!     vrbos => verbosity level
!     ier <= error code, 0 if and only if success
      implicit none
!     --- arguments ---
      integer ncpts, nsrc(3), ndstc(3), cadd
!     ev is at coarse outer points beyond D2
      real*8 ev(-cadd:ndstc(1)+cadd,
     &     -cadd:ndstc(2)+cadd,
     &     -cadd:ndstc(3)+cadd)
!     cp is on a face of D1
      real*8 cp(0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      real*8 sdcrse(ncpts, 0:nsrc(1), 0:nsrc(2), 0:nsrc(3))
      integer vrbos
      integer ier
!     --- locals ---
      integer d2lo(3), d2hi(3), pt
      integer idst, jdst, kdst, isrc, jsrc, ksrc
      integer idir
      real*8 fac, newtrm

      do idir = 1, 3
         d2lo(idir) = 0
         d2hi(idir) = ndstc(idir)
      enddo

      fac = fourth / PI

!     First zero out the coarse outer points.
      do kdst = -cadd, ndstc(3) + cadd
         do jdst = -cadd, ndstc(2) + cadd
            do idst = -cadd, ndstc(1) + cadd
               if (.not. (
     &              (idst .ge. d2lo(1) .and.
     &              idst .le. d2hi(1)) .and.
     &              (jdst .ge. d2lo(2) .and.
     &              jdst .le. d2hi(2)) .and.
     &              (kdst .ge. d2lo(3) .and.
     &              kdst .le. d2hi(3)))) then

                  ev(idst, jdst, kdst) = zero
               endif
            enddo
         enddo
      enddo

!     Loop over all charges on source face.
      do ksrc = 0, nsrc(3)
         do jsrc = 0, nsrc(2)
            do isrc = 0, nsrc(1)

!     Loop over all coarse outer points.
               pt = 0
               do kdst = -cadd, ndstc(3) + cadd
                  do jdst = -cadd, ndstc(2) + cadd
                     do idst = -cadd, ndstc(1) + cadd

                        if (.not. (
     &                       (idst .ge. d2lo(1) .and.
     &                       idst .le. d2hi(1)) .and.
     &                       (jdst .ge. d2lo(2) .and.
     &                       jdst .le. d2hi(2)) .and.
     &                       (kdst .ge. d2lo(3) .and.
     &                       kdst .le. d2hi(3)))) then

                           pt = pt + 1

                           newtrm = cp(isrc, jsrc, ksrc) *
     &                          sdcrse(pt, isrc, jsrc, ksrc)
                           ev(idst, jdst, kdst) =
     &                          ev(idst, jdst, kdst) + newtrm
                        endif

                     enddo
                  enddo
               enddo
!     End of loop over all coarse outer points.
               if (pt .ne. ncpts) then
                  print *, 'evaldirectrectoutercoarse: ',
     &                 'ran through ', pt,
     &                 ' coarse points instead of ', ncpts
                  ier = 2
                  return
               endif

            enddo
         enddo
      enddo

      ier = 0

      return
      end
