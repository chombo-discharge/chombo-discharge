c
c
c     ==================================================================
      subroutine rptt3(ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,imp,impt,bsasdq,
     &                  cmbsasdq,cpbsasdq)
c     ==================================================================
c
c     # Double transverse Riemann solver for the acoustics equations
c     # with varying material properties.

c     # auxN(i,1,:) holds impedance Z
c     # auxN(i,2,:) holds sound speed c
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c
c     #    bsasdq is an array of flux differences that result from a
c     #         transverse splitting (a previous call to rpt3).
c     #         This stands for B^* A^* \Dq but could represent any of
c     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
c     #         and icoor (see below).
c     #         Moreover, each * represents either + or -, as specified by
c     #         imp and impt.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z
c     #      ixyz=2:        y                  z                  x
c     #      ixyz=3:        z                  x                  y
c
c     #    icoor indicates direction in which the transverse solve should
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be
c     #                        split in z into
c     #                           cmbsasdq = C^-B^*A^*\Dq,
c     #                           cpbsasdq = C^+B^*A^*\Dq.
c     #
c     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
c     #                        split in x into
c     #                           cmbsasdq = A^-C^*B^*\Dq,
c     #                           cpbsasdq = A^+C^*B^*\Dq.
c
c     #    The parameters imp and impt are generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients:

c     #    imp =  1 if bsasdq = B^*A^- \Dq, a left-going flux difference
c     #           2 if bsasdq = B^*A^+ \Dq, a right-going flux difference
c     #    impt = 1 if bsasdq = B^-A^* \Dq, a down-going flux difference
c     #           2 if bsasdq = B^+A^* \Dq, an up-going flux difference
c
c     #    aux2(:,:,2) is a 1d slice of the aux array along the row
c     #                 where the data ql, qr lie.
c     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
c     #                 y-like direction
c     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
c     #                z-like direction
c
c
      implicit none

      integer ixyz, icoor, maxm, meqn, mwaves, mbc,mx, maux, imp, impt

      double precision       ql(1-mbc:maxm+mbc, meqn)
      double precision       qr(1-mbc:maxm+mbc, meqn)
      double precision   bsasdq(1-mbc:maxm+mbc, meqn)
      double precision cmbsasdq(1-mbc:maxm+mbc, meqn)
      double precision cpbsasdq(1-mbc:maxm+mbc, meqn)
      double precision     aux1(1-mbc:maxm+mbc, maux, 3)
      double precision     aux2(1-mbc:maxm+mbc, maux, 3)
      double precision     aux3(1-mbc:maxm+mbc, maux, 3)


      integer iuvw, i, i1
      double precision zm, zz, zp, cm, c, cp, a1, a2


c     # set iuvw = 2,3,4, depending on which component of q represents
c     # velocity in the transverse direction in which splitting is to
c     # be performed:
      iuvw = ixyz + icoor
      if (iuvw .gt. 4) then
         iuvw = iuvw-3
      endif

c
      do i = 2-mbc,mx+mbc

c        # The flux difference bsasdq is split into downward moving part
c        # traveling at speed -c relative to the medium below and
c        # an upward moving part traveling at speed +c relative to
c        # the medium above, in the double-transverse direction specified
c        # by icoor.
c
c        # Note that the sum of these parts does not give all of bsasdq
c        # since there is also reflection at the interfaces which decreases
c        # the flux.
c
c
c        # set impendance and sound speed in each row of cells,
c        # by selecting appropriate values from aux arrays depending on
c        # values of icoor, imp, impt:
c
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
c
         if (impt .eq. 1) then
c           # bsasdq propagates to "left" in transverse direction used in rpt3
c           # so we either use auxN(:,:,1) or aux1(:,:,N) for N=1,2,3
c           # depending on icoor:

            if (icoor .eq. 2) then
c           # new double-transverse direction is y-like direction
               zm = aux1(i1,1,1)
               zz = aux2(i1,1,1)
               zp = aux3(i1,1,1)
               cm = aux1(i1,2,1)
               c =  aux2(i1,2,1)
               cp = aux3(i1,2,1)
            else !! (icoor .eq. 3)
c           # new double-transverse direction is z-like direction
               zm = aux1(i1,1,1)
               zz = aux1(i1,1,2)
               zp = aux1(i1,1,3)
               cm = aux1(i1,2,1)
               c =  aux1(i1,2,2)
               cp = aux1(i1,2,3)
            endif
         else
c           # bsasdq propagates to "right" in transverse direction used in rpt3
c           # so we either use auxN(:,:,3) or aux3(:,:,N) for N=1,2,3
c           # depending on icoor:
            if (icoor .eq. 2) then
c           # new double-transverse direction is y-like direction
               zm = aux1(i1,1,3)
               zz = aux2(i1,1,3)
               zp = aux3(i1,1,3)
               cm = aux1(i1,2,3)
               c =  aux2(i1,2,3)
               cp = aux3(i1,2,3)
            else !! (icoor .eq. 3)
c           # new double-transverse direction is z-like direction
               zm = aux3(i1,1,1)
               zz = aux3(i1,1,2)
               zp = aux3(i1,1,3)
               cm = aux3(i1,2,1)
               c =  aux3(i1,2,2)
               cp = aux3(i1,2,3)
            endif
         endif
c

c        # transmitted part of down-going wave:
         a1 = (-bsasdq(i,1) + bsasdq(i,iuvw)*zz) /
     &         (zm + zz)

c        # transmitted part of up-going wave:
         a2 = (bsasdq(i,1) + bsasdq(i,iuvw)*zz) /
     &         (zz + zp)
c
c        # The down-going flux difference cmbsasdq is the product  -c * wave
c
         cmbsasdq(i,1) = cm * a1*zm
         cmbsasdq(i,2) = 0.d0
         cmbsasdq(i,3) = 0.d0
         cmbsasdq(i,4) = 0.d0
         cmbsasdq(i,iuvw) = -cm * a1
c
c        # The up-going flux difference cpbsasdq is the product  c * wave
c
         cpbsasdq(i,1) = cp * a2*zp
         cpbsasdq(i,2) = 0.d0
         cpbsasdq(i,3) = 0.d0
         cpbsasdq(i,4) = 0.d0
         cpbsasdq(i,iuvw) = cp * a2
      enddo

      end
