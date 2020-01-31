C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "CONSTANTS.H"
#include "REAL.H"

      real_t function stretchjacobian(xmap, xlo, xlen, bxwid)
!     Computes Jacobian determinant from Cartesian coordinates.
!
!     => inputs,  <= outputs
!     xmap => cartesian coordinate
!     xlo => lowest of this coordinate
!     xlen => length of this coordinate
!     bxwid => physical box width
      implicit none
!     --- arguments ---
      real_t xmap, xlo, xlen, bxwid
!     --- locals ---
      real_t offctr

!     offctr in range [-1/2:+1/2]
      offctr = (xmap - xlo) / xlen - half
      stretchjacobian = (bxwid*Pi*fourth) / (cos(offctr*Pi*half)**2)

      return
      end
