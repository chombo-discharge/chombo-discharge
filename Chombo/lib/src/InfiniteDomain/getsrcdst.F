      subroutine getsrcdst(srcdst, iev, icp,
     &     n1, s1, s2, ptchsz, npchs, n2c, deg, layer, h,
     &     cpind, cfac, cpow, nterms, vrbos, ier)
c     This function evaluates multipole expansions from patches on a face
c     on the coarse points of another face.
c     ev:  output, evaluation
c     cp:  input, multipole expansions on patches of a face
c     iev:  input, +/- 1,2,3; destination face on which to evaluate
c     icp:  input, +/- 1,2,3; source face with multipole expansions
c     n1:  input, length of domain D1
c     s1:  input, offset from D to D1
c     s2:  input, offset from D1 to D2
c     ptchsz:  input, length of patches
c     npchs:  input, number of patches in each dimension
c     n2c:  input, ceil(n2/ptchsz), number of coarse points in each dimension
c     deg:  input, degree of multipole expansions
c     layer:  input, length of extra layer of coarse points to evaluate
c     h:  input, mesh spacing
c     cpind, cfac, cpow:  input, from subroutine deriv1overr, describing
c                         expressions for partial derivatives of 1/r
c     nterms:  input, number of terms in cpind, cfac, cpow
c     vrbos:  input, verbosity level
c     ier:  output, error code, 0 if and only if success
      implicit none
c     --- arguments ---
      integer iev, icp, n1, s1, s2, ptchsz, npchs, n2c, deg, layer
c     ev is at coarse points of a face of D2, plus a layer
c      real*8 ev(-layer:n2c+layer, -layer:n2c+layer)
c     cp is at centers of patches of a face of D1
c      real*8 cp(0:npchs-1, 0:npchs-1, ((deg+1)*(deg+2))/2)
      real*8 h
      integer nterms
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
      real*8 srcdst(((deg+1)*(deg+2))/2,
     &     -layer:n2c+layer, -layer:n2c+layer, 0:npchs-1, 0:npchs-1)
      integer vrbos
      integer ier
c     --- locals ---
      integer patch1, patch2, dest1, dest2, nc
      integer n2, n
      integer plsrc, pldst, sidsrc, siddst, offsrc, offdst
      integer isrc1, isrc2, idst1, idst2, term, subtrm
      real*8 hbig
      real*8 bsrc(3), bdst(3), src(3), dst(3)
      real*8 xisrc1, xisrc2, xplsrc, r, r2inv, sum
      integer iergpi, dd
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)

c     number of coefficients
      nc = ((deg+1)*(deg+2))/2
c     length of D
      n = n1 - 2*s1
c     length of D2
      n2 = n1 + 2*s2
c     mesh spacing for coarsened level
      hbig = ptchsz * h

c     indlo(term) = min {i: cpind(i) == term}
c     indhi(term) = max {i: cpind(i) == term}
      do term = 1, nterms
         indhi(cpind(term)) = term
      enddo
      do term = nterms, 1, -1
         indlo(cpind(term)) = term
      enddo

c     set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
c     the two perpendicular directions on the source face
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getsrcdst:  illegal icp = ', icp
         ier = 1
         return
      endif

c     set pldst=abs(iev), siddst=sign(iev), and let idst1, idst2 be
c     the two perpendicular directions on the destination face
      call getplaneindices(iev, pldst, siddst, idst1, idst2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getsrcdst:  illegal iev = ', iev
         ier = 2
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc, hilo(siddst), pldst
 110     format (1x, 'source face ', a2, i2, '; dest face ', a2, i2)
      endif

c     Set bsrc(1:3) to be space coordinates of
c     the center of the (0,0) patch on the source face of D1.
c     A face of D1 extends in space
c     from -s1*h to (n+s1)*h = (n1-s1)*h.
c     And the center of the (0, 0) patch is set off from the corner by
c     (ptchsz/2 - offsrc, ptchsz/2 - offsrc) * h.
      offsrc = (npchs * ptchsz - n1) / 2
      bsrc(isrc1) = (ptchsz/2 - offsrc - s1) * h
      bsrc(isrc2) = (ptchsz/2 - offsrc - s1) * h
      if (sidsrc .eq. -1) then
         bsrc(plsrc) = -s1 * h
      elseif (sidsrc .eq. 1) then
         bsrc(plsrc) = (n + s1) * h
      endif
      if (vrbos .gt. 1) then
         print *, 'source is based at ', bsrc(1), bsrc(2), bsrc(3),
     &        ', fixed dimension ', plsrc
      endif

c     src(1:3) holds space coordinates of center of source patch.
c     fixed dimension plsrc; later we'll set
c     src(isrc1) = bsrc(isrc1) + patch1*hbig
c     src(isrc2) = bsrc(isrc2) + patch2*hbig
      src(plsrc) = bsrc(plsrc)

c     Set bdst(1:3) to be space coordinates of
c     the point with evaluation at ev(0, 0).
c     This point is (offdst, offdst) below a corner of D2.
c     Recall that we also evaluate at points with negative
c     indices into ev:  these are points outside D2.
c     A face of D2 extends in space
c     from (-s1-s2)*h to (n+s1+s2)*h = (n1-s1+s2)*h.
      offdst = (n2c * ptchsz - n2) / 2
      bdst(idst1) = (-s1 - s2 - offdst) * h
      bdst(idst2) = (-s1 - s2 - offdst) * h
      if (siddst .eq. -1) then
         bdst(pldst) = (-s1 - s2) * h
      elseif (siddst .eq. 1) then
         bdst(pldst) = (n + s1 + s2) * h
      endif
      if (vrbos .gt. 1) then
         print *, 'dest is based at ', bdst(1), bdst(2), bdst(3),
     &        ', fixed dimension ', pldst
      endif
c     dst(1:3) holds space coordinates of destination point.
c     fixed dimension pldst; later we'll set
c     dst(idst1) = bdst(idst1) + dest1*hbig
c     dst(idst2) = bdst(idst2) + dest2*hbig
      dst(pldst) = bdst(pldst)

c     Loop over all patches on source face.
      do patch2 = 0, npchs-1
         src(isrc2) = bsrc(isrc2) + patch2*hbig
         do patch1 = 0, npchs-1
            src(isrc1) = bsrc(isrc1) + patch1*hbig
c     Patch with index (patch1, patch2) is centered at src(1:3),
c     and has multipole expansion cp(patch1, patch2, 1:nc).
c     print *, 'patch (', patch1, ', ', patch2, ') center ', src
c     Loop over all destination points on destination face.
            do dest2 = -layer, n2c + layer
               dst(idst2) = bdst(idst2) + dest2*hbig
               do dest1 = -layer, n2c + layer
                  dst(idst1) = bdst(idst1) + dest1*hbig
c     Destination point is dst(1:3), and we add to ev(dest1, dest2).
c     xisrc1, xisrc2, xplsrc is displacement vector from src to dst.
c     the multipole expansions are based on powers of xisrc1 and xisrc2.
                  xisrc1 = dst(isrc1) - src(isrc1)
                  xisrc2 = dst(isrc2) - src(isrc2)
                  xplsrc = dst(plsrc) - src(plsrc)
c     r is length of x.
                  r = dsqrt(xisrc1**2 + xisrc2**2 + xplsrc**2)
                  r2inv = 1d0 / (r*r)

c     pow1(k) == xisrc1**k, for k = 0:deg
c     pow2(k) == xisrc2**k, for k = 0:deg
c     powr(k) == r**k, for k = -2*deg-1:0
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
     &                    cfac(subtrm) *
     &                    pow1(cpow(1, subtrm)) *
     &                    pow2(cpow(2, subtrm)) *
     &                    powr(cpow(3, subtrm))
                     enddo
                     srcdst(term, dest1, dest2, patch1, patch2) = sum
                  enddo

c                  do term = 1, nterms
c                     srcdst(cpind(term), dest1, dest2,
c     &                    patch1, patch2) =
c     &                    srcdst(cpind(term), dest1, dest2,
c     &                    patch1, patch2) +
c     &                    cfac(term) *
c     &                    pow1(cpow(1, term)) *
c     &                    pow2(cpow(2, term)) *
c     &                    powr(cpow(3, term))

c                     newtrm =
c     &                    cp(patch1, patch2, cpind(term)) *
c     &                    cfac(term) *
c     &                    (xisrc1**cpow(1, term)) *
c     &                    (xisrc2**cpow(2, term)) *
c     &                    (r**cpow(3, term))

c     Recall that we set sum = ev(dest1, dest2) before adding new terms.
c                  ev(dest1, dest2) = sum

               enddo
            enddo
c     end of loop over destination points
         enddo
      enddo
c     end of loop over patches

      ier = 0

      return
      end



      subroutine getcrsesrcdst(srcdst, icp,
     &     n1, s1, s2, ptchsz, npchs, nbigc, deg, h, cl, refc,
     &     cpind, cfac, cpow, nterms, ncpts, vrbos, ier)
c     This function evaluates multipole expansions from patches on a face
c     of D1 on the coarse points in a region outside of D2.
c     ev:  output, evaluation
c     cp:  input, multipole expansions on patches of a face
c     iev:  input, +/- 1,2,3; destination face on which to evaluate
c     icp:  input, +/- 1,2,3; source face with multipole expansions
c     n1:  input, length of domain D1
c     s1:  input, offset from D to D1
c     s2:  input, offset from D1 to D2
c     ptchsz:  input, length of patches
c     npchs:  input, number of patches in each dimension
c     nbigc:  input, n/refc, number of coarse points in each dimension
c     deg:  input, degree of multipole expansions
c     h:  input, mesh spacing
c     cl:  input, length of extra layer of coarse points to evaluate
c     refc:  input, refinement ratio between fine and coarse
c     cpind, cfac, cpow:  input, from subroutine deriv1overr, describing
c                         expressions for partial derivatives of 1/r
c     nterms:  input, number of terms in cpind, cfac, cpow
c     ncpts:  number of coarse points
c     vrbos:  input, verbosity level
c     ier:  output, error code, 0 if and only if success
      implicit none
c     --- arguments ---
      integer icp, n1, s1, s2, ptchsz, npchs, nbigc, deg, cl, refc
      real*8 h
      integer nterms, ncpts
      integer cpind(nterms)
      integer cpow(3, nterms)
      real*8 cfac(nterms)
c     petermc, 1 mar 2004
c      real*8 srcdst(((deg+1)*(deg+2))/2,
c     &     -cl:nbigc+cl, -cl:nbigc+cl, -cl:nbigc+cl,
c     &     0:npchs-1, 0:npchs-1)
      real*8 srcdst(((deg+1)*(deg+2))/2, ncpts,
     &     0:npchs-1, 0:npchs-1)
      integer vrbos
      integer ier
c     --- locals ---
      integer patch1, patch2, i, j, k, nc, offc
      integer n2, n
      integer plsrc, sidsrc, offsrc
      integer isrc1, isrc2, term, subtrm
      integer d2lo, d2hi
      real*8 hbig, hpatch
      real*8 bsrc(3), src(3), dst(3), x(3)
      real*8 r, r2inv, sum
      integer iergpi, dd, pt
      character*2 hilo(-1:1)
      real*8 pow1(0:deg), pow2(0:deg), powr(0:deg)
      integer indlo(((deg+1)*(deg+2))/2)
      integer indhi(((deg+1)*(deg+2))/2)

c     number of coefficients
      nc = ((deg+1)*(deg+2))/2
c     length of D
      n = n1 - 2*s1
c     length of D2
      n2 = n1 + 2*s2
c     space between patches
      hpatch = ptchsz * h
c     mesh spacing for coarsened level
      hbig = refc * h

c     indhi(term) = max {i: cpind(i) == term}
      do term = 1, nterms
         indhi(cpind(term)) = term
      enddo

c     indlo(term) = min {i: cpind(i) == term}
      do term = nterms, 1, -1
         indlo(cpind(term)) = term
      enddo

c     set plsrc=abs(icp), sidsrc=sign(icp), and let isrc1, isrc2 be
c     the two perpendicular directions on the source face
      call getplaneindices(icp, plsrc, sidsrc, isrc1, isrc2, iergpi)
      if (iergpi .ne. 0) then
         print *, 'getsrcdst:  illegal icp = ', icp
         ier = 1
         return
      endif

      if (vrbos .gt. 0) then
         hilo(-1) = 'LO'
         hilo(1) = 'HI'
         print 110, hilo(sidsrc), plsrc
 110     format (1x, 'source face ', a2, i2, ' onto coarse points')
      endif

c     Set bsrc(1:3) to be space coordinates of
c     the center of the (0,0) patch on the source face of D1.
c     A face of D1 extends in space
c     from -s1*h to (n+s1)*h = (n1-s1)*h.
c     And the center of the (0, 0) patch is set off from the corner by
c     (ptchsz/2 - offsrc, ptchsz/2 - offsrc) * h.
      offsrc = (npchs * ptchsz - n1) / 2
      bsrc(isrc1) = (ptchsz/2 - offsrc - s1) * h
      bsrc(isrc2) = (ptchsz/2 - offsrc - s1) * h
      if (sidsrc .eq. -1) then
         bsrc(plsrc) = -s1 * h
      elseif (sidsrc .eq. 1) then
         bsrc(plsrc) = (n + s1) * h
      endif
      if (vrbos .gt. 1) then
         print *, 'source is based at ', bsrc(1), bsrc(2), bsrc(3),
     &        ', fixed dimension ', plsrc
      endif

c     src(1:3) holds space coordinates of center of source patch.
c     fixed dimension plsrc; later we'll set
c     src(isrc1) = bsrc(isrc1) + patch1*hbig
c     src(isrc2) = bsrc(isrc2) + patch2*hbig
      src(plsrc) = bsrc(plsrc)

c     Recall that D2 has range -s1-s2:n+s1+s2 in each dimension.
c     Also, refc divides n.
      offc = (s1 + s2) / refc
c     d2lo is the lowest coarse coordinate of a coarse point inside D2
      d2lo = -offc
c     d2hi is the highest coarse coordinate of a coarse point inside D2
      d2hi = nbigc + offc

      pt = 0
c     Loop over all destination points
      do k = -cl, nbigc + cl
         dst(3) = k*hbig
         do j = -cl, nbigc + cl
            dst(2) = j*hbig
            do i = -cl, nbigc + cl
c     check that we're outside D2:
c     at least one coordinate must be outside the range d2lo:d2hi
               if (.not. (
     &              (i .ge. d2lo .and. i .le. d2hi) .and.
     &              (j .ge. d2lo .and. j .le. d2hi) .and.
     &              (k .ge. d2lo .and. k .le. d2hi))) then

                  pt = pt + 1

                  dst(1) = i*hbig

c     Loop over all patches on source face.
                  do patch2 = 0, npchs-1
                     src(isrc2) = bsrc(isrc2) + patch2*hpatch
                     do patch1 = 0, npchs-1
                        src(isrc1) = bsrc(isrc1) + patch1*hpatch
c     Patch with index (patch1, patch2) is centered at src(1:3),
c     and has multipole expansion cp(patch1, patch2, 1:nc).
c     print *, 'patch (', patch1, ', ', patch2, ') center ', src

c     Destination point is dst(1:3), and we add to ev(i, j, k).
c     x is displacement vector from src to dst.
c     The multipole expansions are based on powers of x(isrc1) and x(isrc2).
                        x(1) = dst(1) - src(1)
                        x(2) = dst(2) - src(2)
                        x(3) = dst(3) - src(3)
c     r is length of x.
                        r = dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
                        r2inv = 1d0 / (r*r)

c     pow1(k) == x(isrc1)**k, for k = 0:deg
c     pow2(k) == x(isrc2)**k, for k = 0:deg
c     powr(k) == 1/r**(2*k+1), for k = 0:deg
                        pow1(0) = 1d0
                        pow2(0) = 1d0
                        powr(0) = 1d0 / r
                        do dd = 1, deg
                           pow1(dd) = x(isrc1) * pow1(dd-1)
                           pow2(dd) = x(isrc2) * pow2(dd-1)
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

c                           srcdst(term, i, j, k, patch1, patch2) =
c     &                          sum
                           srcdst(term, pt, patch1, patch2) = sum

                        enddo

c     Recall that we set sum = ev(dest1, dest2) before adding new terms.
c     ev(dest1, dest2) = sum
                     enddo
                  enddo
c     end of loop over patches
               endif
c     end of check whether in range d2lo:d2hi
            enddo
         enddo
      enddo
c     end of loop over destination points

      ier = 0

      return
      end



      subroutine getplaneindices(ind, plane, side, i1, i2, ier)
c     input:
c     ind:  +/- 1, 2, 3
c     outputs:
c     plane:  which dimension is constant (1, 2, 3)
c     side:  -1 for low, +1 for high
c     i1, i2:  the two variable dimensions (1, 2, 3 but not plane)
c     ier:  return code, 0 if and only if success
      implicit none
      integer ind, plane, side, i1, i2, ier
      integer iabs

      plane = iabs(ind)
      side = ind / plane

      if (side .ne. 1 .and. side .ne. -1) then
c     illegal side
         ier = 1
         return
      endif

      if (plane .eq. 1) then
c     plane x = const
         i1 = 2
         i2 = 3
      elseif (plane .eq. 2) then
c     plane y = const
         i1 = 1
         i2 = 3
      elseif (plane .eq. 3) then
c     plane z = const
         i1 = 1
         i2 = 2
      else
c     illegal plane
         ier = 2
         return
      endif

      ier = 0

      return
      end
