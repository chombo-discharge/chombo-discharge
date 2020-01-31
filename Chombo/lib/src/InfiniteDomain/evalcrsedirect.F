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
      subroutine evalcrsedirectfaces(ev, q,
     &     n, s1, s2, nbigc, h, cl, refc, vrbos, ier)
      implicit none
!     --- arguments ---
      integer n, s1, s2, nbigc, cl, refc
!     ev is at coarse points of D2, plus a layer
      real*8 ev(-cl:nbigc+cl, -cl:nbigc+cl, -cl:nbigc+cl)
!     q is at points of faces of D1
!     petermc, 8 Jan 2004:  last two dimensions are x/y/z and lo/hi
      real*8 q(-s1:n+s1, -s1:n+s1, 1:3, 0:1)
      real*8 h
      integer vrbos
      integer ier
!     --- locals ---
      integer iside, idim, icp, ierf

      call zerocube(ev, nbigc + 2*cl)
      do iside = 0, 1
         do idim = 1, 3
            icp = (2*iside - 1) * idim
            call evalcrsedirect(ev, q(-s1, -s1, idim, iside), icp,
     &           n, s1, s2, nbigc, h, cl, refc, vrbos, ierf)
            if (ierf .ne. 0) then
               ier = 200
               if (iside .eq. 0) then
                  ier = ier + (10 - idim)
               else
                  ier = ier + idim
               endif
               return
            endif
         enddo
      enddo

      ier = 0

      return
      end


      subroutine evalcrsedirect(ev, q, icp,
     &     n, s1, s2, nbigc, h, cl, refc, vrbos, ier)
      implicit none
!     --- arguments ---
      integer icp, n, s1, s2, nbigc, cl, refc
!     ev is at coarse points of D2, plus a layer
      real*8 ev(-cl:nbigc+cl, -cl:nbigc+cl, -cl:nbigc+cl)
!     q is at points of a face of D1
      real*8 q(-s1:n+s1, -s1:n+s1)
      real*8 h
      integer vrbos, ier
!     --- locals ---
      integer plsrc, sidsrc, isrc1, isrc2, iergpi, iergi, n1
      integer i, j, k, src1, src2, d2lo, d2hi, offc, intrul
      real*8 sum, xsrc(3), xdst(3), r, pi, hbig, wt(-s1:n+s1)
      real*8 bsrc(3), fac, newtrm, xdif(3)
      character*2 hilo(-1:1)
      real*8 dsqrt, datan

      n1 = n + 2*s1
!     mesh spacing for coarsened level
      hbig = refc * h

!     Set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
!     the two perpendicular directions on the source face.
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         if (vrbos .gt. 0) then
            print *, 'evalcrsempc:  illegal icp = ', icp
         endif
         ier = 1
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 130, hilo(sidsrc), plsrc
 130     format (1x, 'source face ', a2, i2, ' onto coarse points')
      endif

!     Set bsrc(1:3) to be space coordinates of
!     the point indexed by (0, 0) on the source face of D1,
!     where indices are in the range (-s1:n+s1, -s1:n+s1).
!     In the dimensions isrc1 and isrc2, a face of D1 extends
!     in space from -s1*h to (n+s1)*h = (n1-s1)*h.
!     Space coordinates of point indexed by (isrc1, isrc2) are
!     xsrc(isrc1) = bsrc(isrc1) + isrc1 * h
!     xsrc(isrc2) = bsrc(isrc2) + isrc2 * h
!     xsrc(plsrc) = bsrc(plsrc)
      bsrc(isrc1) = 0d0
      bsrc(isrc2) = 0d0
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

      pi = 4 * datan(1d0)
      fac = 0.25d0 / pi

!     1 trapezoidal, 2 Simpson, 4 Bode, 8 Newton-Cotes-nine
      intrul = 4
!     Bode's rule here (intrul=4) is an improvement on
!     Simpson's (intrul=2) when n is big
!     -- can get 1/30 of minimum for n=512 and far faces.
!     Changing to intrul = 8 doesn't buy you anything.
      call getintegral1dweights(wt, -s1, n+s1, h, intrul, iergi)

!     xsrc(1:3) contains space coordinates of source point.
!     xsrc(plsrc) is fixed
      xsrc(plsrc) = bsrc(plsrc)

!     Recall that D2 has range -s1-s2:n+s1+s2 in each dimension.
!     Also, refc divides n.
      offc = (s1 + s2) / refc
!     d2lo is the lowest coarse coordinate of a coarse point inside D2
      d2lo = -offc
!     d2hi is the highest coarse coordinate of a coarse point inside D2
      d2hi = nbigc + offc

!     Loop over all destination points
!     xdst(1:3) contains space coordinates of destination point.
      do k = -cl, nbigc + cl
         xdst(3) = k*hbig
         do j = -cl, nbigc + cl
            xdst(2) = j*hbig
            do i = -cl, nbigc + cl
!     check that we're outside D2
               if (.not.
     &              (i .ge. d2lo .and. i .le. d2hi) .and.
     &              (j .ge. d2lo .and. j .le. d2hi) .and.
     &              (k .ge. d2lo .and. k .le. d2hi)) then

                  xdst(1) = i*hbig

                  sum = ev(i, j, k)

!     xdif(1:3) is displacement vector from xsrc(1:3) to xdst(1:3).
                  xdif(plsrc) = xdst(plsrc) - xsrc(plsrc)

!     Loop over all points on source face.
!     Source point is xsrc(1:3).
                  do src2 = -s1, n+s1
                     xsrc(isrc2) = bsrc(isrc2) + src2 * h
                     xdif(isrc2) = xdst(isrc2) - xsrc(isrc2)

                     do src1 = -s1, n+s1
                        xsrc(isrc1) = bsrc(isrc1) + src1 * h
                        xdif(isrc1) = xdst(isrc1) - xsrc(isrc1)

!     r is length of xdif(1:3).
                        r = dsqrt(xdif(1) * xdif(1)
     &                       + xdif(2) * xdif(2)
     &                       + xdif(3) * xdif(3))

                        newtrm = fac * wt(src1) * wt(src2) *
     &                       q(src1, src2) / r

                        sum = sum + newtrm
                     enddo
                  enddo
!     Recall that we had set sum = ev(i, j, k) before adding new terms.
                  ev(i, j, k) = sum
!     end of loop over source points
               endif
            enddo
         enddo
      enddo
!     end of loop over destination points

      ier = 0

      return
      end
