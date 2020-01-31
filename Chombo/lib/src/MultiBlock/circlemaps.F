C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "CONSTANTS.H"
#include "REAL.H"

      real_t function sectoroffset(iblock)
!     Returns angular offset of block.
!
!     block 1:   -Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4]
!     block 2:    Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + Pi/2
!     block 3:  3*Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + 2*Pi/2
!     block 4:  5*Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + 3*Pi/2
      implicit none
!     --- arguments ---
      integer iblock
!     --- locals ---

      sectoroffset = (iblock - 1) * (Pi*half)

      return
      end


      subroutine spokesmap(xphys, yphys, angmap, radmap,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     xphys <= physical x coordinate
!     yphys <= physical y coordinate
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t xphys, yphys
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
!     --- locals ---
      real_t amoff, radoff, r0, r1, theta

!     amoff in interval [-1/2:+1/2], with amoff=-1/2 at angmap=anglo
      amoff = (angmap - anglo) / anglen - half
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      call getrthetaside(r0, theta, amoff, bxwid)
      theta = theta + thoff
      r1 = outr
      call getrec(xphys, yphys, r0, r1, radoff, theta, xorg, yorg)

      return
      end


      subroutine equalanglesmap(xphys, yphys, angmap, radmap,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     xphys <= physical x coordinate
!     yphys <= physical y coordinate
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t xphys, yphys
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
!     --- locals ---
      real_t amoff, radoff, r0, r1, theta, th

!     amoff in interval [-1/2:+1/2], with amoff=-1/2 at angmap=anglo
      amoff = (angmap - anglo) / anglen - half
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      th = amoff * (Pi*half)
      r0 = bxwid * half / cos(th)
      r1 = outr
      theta = th + thoff

      call getrec(xphys, yphys, r0, r1, radoff, theta, xorg, yorg)

      return
      end


      subroutine transitionmap(xphys, yphys, angmap, radmap,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     xphys <= physical x coordinate
!     yphys <= physical y coordinate
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t xphys, yphys
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
!     --- locals ---
      real_t amoff, radoff, r0, r1, theta, th0, th1

!     amoff in interval [-1/2:+1/2], with amoff=-1/2 at angmap=anglo
      amoff = (angmap - anglo) / anglen - half
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      call getrthetaside(r0, th0, amoff, bxwid)
      r1 = outr
      th1 = amoff * (Pi*half)
      theta = (one - radoff)*th0 + radoff*th1 + thoff

      call getrec(xphys, yphys, r0, r1, radoff, theta, xorg, yorg)

      return
      end


      subroutine getrthetaside(r, theta, amoff, bxwid)
!     Gets polar coordinates.
!
!     => inputs,  <= outputs
!     r <= radius
!     theta <= angle
!     amoff => angular variable offset
!     bxwid => physical box width
      implicit none
!     --- arguments ---
      real_t r, theta, amoff, bxwid
!     --- locals ---
      real_t hyp

      hyp = sqrt(amoff*amoff + fourth)
      theta = asin(amoff / hyp)
      r = bxwid * hyp

      return
      end


      subroutine getrec(xphys, yphys, r0, r1, radmap, theta,
     &     xorg, yorg)
      implicit none
!     --- arguments ---
      real_t xphys, yphys, r0, r1, radmap, theta, xorg, yorg
!     --- locals ---
      real_t r

      r = (one - radmap)*r0 + radmap*r1
      xphys = xorg + r*cos(theta)
      yphys = yorg + r*sin(theta)

      return
      end


!---- inverse maps


      subroutine spokesinvmap(angmap, radmap, xphys, yphys,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr, angsgn, radsgn)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     angmap <= cartesian angular coordinate
!     radmap <= cartesian radial coordinate
!     xphys => physical x coordinate
!     yphys => physical y coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
!     angsgn => sign by which to multiply angmap
!     radsgn => sign by which to multiply radmap
      implicit none
!     --- arguments ---
      real_t angmap, radmap
      real_t xphys, yphys, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
      integer angsgn, radsgn
!     --- locals ---
      real_t r, theta, xoff, rooter, r0, amoff, radoff

      call getpolar(r, theta, xphys, yphys, xorg, yorg, thoff)

!     theta in range [-Pi/4:+Pi/4].
!     Solve for xoff, in range [-1/2:+1/2].
      xoff = tan(theta) * half
      rooter = sqrt(xoff * xoff + fourth)
      r0 = bxwid * rooter

!     amoff, radoff in range [0:1].
      amoff = half + xoff
      radoff = (r - r0) / (outr - r0)

      angmap = angsgn * (anglo + anglen * amoff)
      radmap = radsgn * (radlo + radlen * radoff)

      return
      end


      subroutine equalanglesinvmap(angmap, radmap, xphys, yphys,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr, angsgn, radsgn)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     angmap <= cartesian angular coordinate
!     radmap <= cartesian radial coordinate
!     xphys => physical x coordinate
!     yphys => physical y coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
!     angsgn => sign by which to multiply angmap
!     radsgn => sign by which to multiply radmap
      implicit none
!     --- arguments ---
      real_t angmap, radmap
      real_t xphys, yphys, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
      integer angsgn, radsgn
!     --- locals ---
      real_t r, theta, r0, amoff, radoff

      call getpolar(r, theta, xphys, yphys, xorg, yorg, thoff)

!     theta in range [-Pi/4:+Pi/4].
      r0 = bxwid * half / cos(theta)

!     amoff, radoff in range [0:1].
      amoff = half + theta / (Pi * half)
      radoff = (r - r0) / (outr - r0)

      angmap = angsgn * (anglo + anglen * amoff)
      radmap = radsgn * (radlo + radlen * radoff)

      return
      end


      subroutine transitioninvmap(angmap, radmap, xphys, yphys,
     &     anglo, radlo, anglen, radlen, xorg, yorg,
     &     thoff, bxwid, outr, angsgn, radsgn)
!     Computes physical coordinates from mapped coordinates.
!
!     => inputs,  <= outputs
!     angmap <= cartesian angular coordinate
!     radmap <= cartesian radial coordinate
!     xphys => physical x coordinate
!     yphys => physical y coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     xorg => physical x coordinate of origin of circle
!     yorg => physical y coordinate of origin of circle
!     thoff => angular offset of block
!     bxwid => physical box width
!     outr => radius of circle
!     angsgn => sign by which to multiply angmap
!     radsgn => sign by which to multiply radmap
      implicit none
!     --- arguments ---
      real_t angmap, radmap
      real_t xphys, yphys, anglo, radlo, anglen, radlen
      real_t xorg, yorg, thoff, bxwid, outr
      integer angsgn, radsgn
!     --- locals ---
      real_t r, theta, xoff, rooter, r0, amoff, radoff
      real_t xoffps, xfun, xfund, deltax, thresh, eps
      integer itn, maxitn
      parameter (maxitn=10)
#ifdef CH_USE_DOUBLE
      real_t dlamch
      external dlamch
#else
      real_t slamch
      external slamch
#endif

      call getpolar(r, theta, xphys, yphys, xorg, yorg, thoff)

!     Solve for xoff, in range -1/2:1/2.
!     petermc, 19 May 2008, changed from thoff to theta.
      xoff = theta / (Pi * half)
#ifdef CH_USE_DOUBLE
      eps = dlamch('Epsilon')
#else
      eps = slamch('Epsilon')
#endif
!     if this is 2*eps then I get hang
      thresh = 4*eps
      itn = 0
!     Newton's method
      do while (.true.)
         xoffps = one + four * xoff*xoff
         rooter = sqrt(xoffps)
         xfun = outr*theta - (half*Pi)*r*xoff
     &        + bxwid*fourth*(-two*theta + Pi*xoff)*rooter
     &        + (r - outr)*asin(two*xoff/rooter)
         xfund = (-two*rooter * (bxwid*theta*xoff
     &        + (outr - r)/rooter)
     &        + Pi*fourth*(bxwid*rooter*
     &        (one + eight*xoff*xoff) - two*r*xoffps))
     &        / xoffps
         deltax = xfun / xfund
         if (abs(deltax) .lt. thresh) exit
         itn = itn + 1
         xoff = xoff - deltax
         if (itn .ge. maxitn) then
            print *, 'TRANSITIONINVMAP warning:  hang'
            exit
         endif
      enddo
!     xoff in range [-1/2 : 1/2];
!     hence, rooter in range [1/2 : 1/sqrt(2)].
      xoffps = xoff*xoff + fourth
      rooter = sqrt(xoffps)
!     r0 in range [bxwid/2 : bxwid/sqrt(2)]
      r0 = bxwid * rooter
!     xoff in range [-1/2 : 1/2];
!     r in range [r0 : outr].
!     Hence, amoff, radoff in range [0 : 1].
      amoff = half + xoff
      radoff = (r - r0) / (outr - r0)

      angmap = angsgn * (anglo + anglen * amoff)
      radmap = radsgn * (radlo + radlen * radoff)

      return
      end


      subroutine getpolar(r, theta, xphys, yphys, xorg, yorg, thoff)
      implicit none
!     --- arguments ---
      real_t r, theta, xphys, yphys, xorg, yorg, thoff
!     --- locals ---
      real_t xtrans, ytrans

      xtrans = xphys - xorg
      ytrans = yphys - yorg

      r = sqrt(xtrans*xtrans + ytrans*ytrans)

!     block 1:   -Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4]
!     block 2:    Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + Pi/2
!     block 3:  3*Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + 2*Pi/2
!     block 4:  5*Pi/4 + [0:Pi/2] = [-Pi/4:+Pi/4] + 3*Pi/2
!     Put theta in range 0:Pi/2 within its block.
      theta = atan2(ytrans, xtrans)
!     theta = th - (2*iblock - 3) * fourth * Pi
      theta = theta - thoff + fourth * Pi
!     Put theta in range 0:2*Pi.
      theta = mod(theta + two * Pi, two * Pi)
!     Put theta in range -Pi/4:+Pi/4.
      theta = theta - fourth * Pi

      return
      end


!---- Jacobians


      real_t function spokesjacobian(angmap, radmap,
     &     anglo, radlo, anglen, radlen, bxwid, outr)
!     Computes Jacobian determinant from Cartesian coordinates.
!
!     => inputs,  <= outputs
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t bxwid, outr
!     --- locals ---
      real_t amoff, radoff, parlen

!     amoff in interval [0:1], with amoff=0 at angmap=anglo
      amoff = (angmap - anglo) / anglen
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      parlen = one + two*amoff*(amoff - one)
      spokesjacobian =
     &     - ((one - radoff)*(bxwid**2)*parlen - two*radoff*(outr**2)
     &     + (two*radoff - one)*bxwid*outr * sqrt(two*parlen))
     &     / (two * parlen)

      return
      end


      real_t function equalanglesjacobian(angmap, radmap,
     &     anglo, radlo, anglen, radlen, bxwid, outr)
!     Computes Jacobian determinant from Cartesian coordinates.
!
!     => inputs,  <= outputs
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t bxwid, outr
!     --- locals ---
      real_t amoff, radoff, ccos

!     amoff in interval [0:1], with amoff=0 at angmap=anglo
      amoff = (angmap - anglo) / anglen
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      ccos = cos((two*amoff - one)*Pi*fourth)
      equalanglesjacobian =
     &     - Pi * (bxwid - two*outr*ccos) *
     &     (bxwid + radoff*(two*outr*ccos - bxwid)) /
     &     (eight*ccos*ccos)

      return
      end


      real_t function transitionjacobian(angmap, radmap,
     &     anglo, radlo, anglen, radlen, bxwid, outr)
!     Computes Jacobian determinant from Cartesian coordinates.
!
!     => inputs,  <= outputs
!     angmap => cartesian angular coordinate
!     radmap => cartesian radial coordinate
!     anglo => lowest angular coordinate in sector
!     radlo => lowest radial coordinate in sector
!     anglen => angular coordinate length of sector
!     radlen => radial coordinate length of sector
!     bxwid => physical box width
!     outr => radius of circle
      implicit none
!     --- arguments ---
      real_t angmap, radmap, anglo, radlo, anglen, radlen
      real_t bxwid, outr
!     --- locals ---
      real_t amoff, radoff, parlen, sqpl, sq2pl, radom1, amo2m1

!     amoff in interval [0:1], with amoff=0 at angmap=anglo
      amoff = (angmap - anglo) / anglen
!     radoff in interval [0:1], with radoff=0 at radmap=radlo
      radoff = (radmap - radlo) / radlen

      parlen = one + two*amoff*(amoff - one)
      sqpl = sqrt(parlen)
      sq2pl = sqrt(two*parlen)

      radom1 = radoff - one
      amo2m1 = two*amoff - one

      transitionjacobian =
     &     ( (bxwid*sq2pl*radom1 - two*outr*radoff)
     &     * (four*outr*(two*radom1 - Pi*parlen*radoff)
     &     + sqrt(two)*bxwid*
     &     ((four*Pi*sqpl - eight*radom1/sqpl)*amoff*(amoff - one)
     &     - four*radom1/sqpl + Pi*(radoff + one)*sqpl) +
     &     four*bxwid*amo2m1*radom1*sq2pl
     &     *asin(amo2m1/sq2pl))) / (sixteen*parlen)

      return
      end


