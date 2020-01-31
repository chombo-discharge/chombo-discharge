C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "REAL.H"
#include "CONSTANTS.H"

      real_t function BIzonalvelocity(phi)
      implicit none
#include "SWEosCommon.fh"
#include "SWEBICommon.fh"
!     --- arguments ---
      real_t phi
!     --- locals ---

      BIzonalvelocity = zero
      if ((phi0 .lt. phi) .and. (phi .lt. phi1)) then
         BIzonalvelocity =
     &        (umax/en) * exp(one/((phi - phi0)*(phi - phi1)))
      endif

      return
      end


      real_t function BIintegrand(phi)
      implicit none
#include "SWEosCommon.fh"
#include "SWEBICommon.fh"
!     --- arguments ---
      real_t phi
!     --- locals ---
      real_t f, u, BIzonalvelocity

      f = two * omega * sin(phi)
      
!     BIzonalvelocity(phi) is zero outside [phi0, phi1],
      BIintegrand = zero
      if ((phi0 .lt. phi) .and. (phi .lt. phi1)) then
         u = BIzonalvelocity(phi)
         BIintegrand = u * (f + tan(phi) * u)
      endif

      return
      end


      real_t function BIintegrand2(phi)
      implicit none
#include "SWEosCommon.fh"
#include "SWEBICommon.fh"
!     --- arguments ---
      real_t phi
!     --- locals ---
      real_t rombergintegral, romint
      external BIintegrand

      romint = rombergintegral(BIintegrand, -half*Pi, phi, 14)
      BIintegrand2 = half * cos(phi) * romint

      return
      end


      real_t function rombergintegral(f, a, b, itns)
!
!     => inputs
!     f       => integrand
!     a       => lower limit of integral
!     b       => upper limit of integral
!     itns    => number of iterations; try 10
!
      implicit none
!     --- arguments ---
      real_t f, a, b
      integer itns
!     --- locals ---
      real_t rombergintegralinner
      integer nintvl

      nintvl = 2**(itns-1)
      rombergintegral = rombergintegralinner(f, a, b, itns, nintvl)

      return
      end


      real_t function rombergintegralinner(f, a, b, itns, nintvl)
      implicit none
!     --- arguments ---
      real_t f, a, b
      integer itns, nintvl
!     --- locals ---
      real_t romall(0:nintvl), h, pt, stepsz
      real_t rom1(itns), rom2(itns), romsum
      integer st, i, k

      h = b - a
      stepsz = h / (nintvl * one)
      do k = 0, nintvl
         pt = a + k * stepsz
         romall(k) = f(pt)
      enddo
      rom1(1) = h*(romall(0) + romall(nintvl))*half
      st = nintvl
      do i = 2, itns
         romsum = zero
         do k = st/2, nintvl, st
            romsum = romsum + romall(k)
         enddo
         rom2(1) = (rom1(1) + h*romsum)*half
         do k = 1, i-1
            rom2(k+1) = ((four**k)*rom2(k) - rom1(k)) /
     &           (four**k - one)
         enddo
         do k = 1, i
            rom1(k) = rom2(k)
         enddo
         h = half * h
         st = st / 2
      enddo
      rombergintegralinner = rom1(itns)
      
      return
      end
