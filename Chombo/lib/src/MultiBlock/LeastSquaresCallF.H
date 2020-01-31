#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEASTSQUARESCALLF_H_
#define _LEASTSQUARESCALLF_H_

#include "FORT_PROTO.H"
#include "REAL.H"

extern "C"
{
#define FORT_LSWORKSIZE FORTRAN_NAME(LSWORKSIZE, lsworksize)
  int FORT_LSWORKSIZE
  (
   const int* neqns,
   const int* nvars
   );

#define FORT_SOLVELS FORTRAN_NAME(SOLVELS, solvels)
  void FORT_SOLVELS
  (
   const int* neqns,
   const int* nvars,
   Real* mat,
   Real* vec,
   const int* lwork,
   Real* work
   );
}

#endif
