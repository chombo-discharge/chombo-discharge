C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "CONSTANTS.H"
#include "REAL.H"

      subroutine spheremap(xphys, yphys, zphys, rmap, amap, bmap,
     &     rlo, alo, blo, rlen, alen, blen, xorg, yorg, zorg,
     &     iblock, bxwid, outr)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     xphys <= physical x coordinate
!     yphys <= physical y coordinate
!     zphys <= physical z coordinate
!     rmap => mapped RADIAL coordinate
!     amap => mapped ADIR coordinate
!     bmap => mapped BDIR coordinate
!     rlo => lowest in RADIAL coordinate
!     alo => lowest in ADIR coordinate
!     blo => lowest in BDIR coordinate
!     rlen => length in RADIAL coordinate
!     alen => length in ADIR coordinate
!     blen => length in BDIR coordinate
!     xorg => x component of physical x
!     yorg => y component of physical y
!     zorg => z component of physical z
!     iblock => index of block, in 0:5
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      integer iblock
      real_t xphys, yphys, zphys, rmap, amap, bmap
      real_t rlo, alo, blo, rlen, alen, blen
      real_t xorg, yorg, zorg, bxwid, outr
!     --- locals ---
      real_t rrat, arat, brat
      real_t a, b, r, ee, rff, rats(0:2)
      integer vars(0:2, 0:5), signs(0:2, 0:5)
      data vars / 0,2,1, 2,0,1, 0,2,1, 2,0,1, 1,2,0, 1,2,0 /
      data signs / 1,1,1, -1,1,1, -1,-1,1, 1,-1,1, -1,1,1, 1,1,-1 /

      rrat = (rmap - rlo) / rlen
      arat = (amap - alo) / alen
      brat = (bmap - blo) / blen

      r = (one - rrat)*bxwid*half + rrat*outr
      a = two*arat - one
      b = two*brat - one

      ee = one + a*a + b*b
      rff = r / sqrt(((outr - r) + (r - bxwid*half)*ee) /
     &     (outr - bxwid*half))

      rats(0) = one
      rats(1) = a
      rats(2) = b

      xphys = xorg + signs(0, iblock) * rats(vars(0, iblock)) * rff
      yphys = yorg + signs(1, iblock) * rats(vars(1, iblock)) * rff
      zphys = zorg + signs(2, iblock) * rats(vars(2, iblock)) * rff

      return
      end


      subroutine sphereinvmap(rmap, amap, bmap, xphys, yphys, zphys,
     &     rlo, alo, blo, rlen, alen, blen, xorg, yorg, zorg,
     &     iblock, bxwid, outr, rsign, asign, bsign)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     rmap <= mapped RADIAL coordinate
!     amap <= mapped ADIR coordinate
!     bmap <= mapped BDIR coordinate
!     xphys => physical x coordinate
!     yphys => physical y coordinate
!     zphys => physical z coordinate
!     rlo => lowest in RADIAL coordinate
!     alo => lowest in ADIR coordinate
!     blo => lowest in BDIR coordinate
!     rlen => length in RADIAL coordinate
!     alen => length in ADIR coordinate
!     blen => length in BDIR coordinate
!     xorg => x component of physical x
!     yorg => y component of physical y
!     zorg => z component of physical z
!     iblock => index of block, in 0:5
!     bxwid => physical box width
!     outr => radius of circle
!     rsign => sign of RADIAL
!     asign => sign of ADIR
!     bsign => sign of BDIR
      implicit none
!     --- arguments ---
      integer iblock
      real_t rmap, amap, bmap, xphys, yphys, zphys
      real_t rlo, alo, blo, rlen, alen, blen
      real_t xorg, yorg, zorg, bxwid, outr
      integer rsign, asign, bsign
!     --- locals ---
      real_t xtrans, ytrans, ztrans, a, b, ee
      real_t rrat, arat, brat
      real_t r2, aq, bq, cq, r, rats(0:2)
      integer numa(0:5), numb(0:5), denom(0:5), signa(0:5), signb(0:5)
      data numa / 2, 2, 2, 2, 0, 0 /
      data numb / 1, 0, 1, 0, 1, 1 /
      data denom / 0, 1, 0, 1, 2, 2 /
      data signa / 1, 1, -1, -1, -1, -1 /
      data signb / 1, -1, 1, -1, 1, -1 /

      xtrans = xphys - xorg
      ytrans = yphys - yorg
      ztrans = zphys - zorg
      r2 = xtrans*xtrans + ytrans*ytrans + ztrans*ztrans

      rats(0) = xtrans
      rats(1) = ytrans
      rats(2) = ztrans

      a = signa(iblock) * rats(numa(iblock)) / rats(denom(iblock))
      b = signb(iblock) * rats(numb(iblock)) / rats(denom(iblock))

!     coefficients of quadratic
      ee = one + a*a + b*b
      aq = ee * (outr - bxwid * half)
      bq = - r2 * (a*a + b*b)
      cq = r2 * (bxwid*half*ee - outr)
!     Take only the positive quadratic root, because aq*cq < 0,
!     and hence bq^2 - 4*aq*cq > bq^2.
      r = (-bq + sqrt(bq*bq - four*aq*cq)) / (two*aq)

      rrat = (r - bxwid*half) / (outr - bxwid*half)
      arat = (one + a) * half
      brat = (one + b) * half

      rmap = rsign * (rlo + rrat * rlen)
      amap = asign * (alo + arat * alen)
      bmap = bsign * (blo + brat * blen)

      return
      end


      real_t function spherejacobian(rmap, amap, bmap,
     &     rlo, alo, blo, rlen, alen, blen,
     &     bxwid, outr)
!     Computes Jacobian determinant from mapped coordinates.
!
!     => inputs,  <= outputs
!     rmap => mapped RADIAL coordinate
!     amap => mapped ADIR coordinate
!     bmap => mapped BDIR coordinate
!     rlo => lowest in RADIAL coordinate
!     alo => lowest in ADIR coordinate
!     blo => lowest in BDIR coordinate
!     rlen => length in RADIAL coordinate
!     alen => length in ADIR coordinate
!     blen => length in BDIR coordinate
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t rmap, amap, bmap, bxwid, outr
      real_t rlo, alo, blo, rlen, alen, blen
!     --- locals ---
      real_t rrat, arat, brat
      real_t xyp

      rrat = (rmap - rlo) / rlen
      arat = (amap - alo) / alen
      brat = (bmap - blo) / blen

      xyp = rrat*(rrat - one) + arat*(arat - one)
!     5 Mar 2008:  negated
      spherejacobian =
     &     - ((bxwid * (brat - one) - two*outr*brat)**2 *
     &     (two*(bxwid*(one + xyp) - outr)
     &     + (bxwid - two*outr)*(one + two*xyp)*brat)) /
     &     (two * sqrt(one + (two + four*xyp)*brat)**5)

      return
      end
