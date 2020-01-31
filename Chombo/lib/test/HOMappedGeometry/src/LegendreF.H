#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEGENDREF_H_
#define _LEGENDREF_H_


#include "FORT_PROTO.H"
#include "REAL.H"

extern "C"
{
#define FORT_LEGENDRECOEFFS FORTRAN_NAME(LEGENDRECOEFFS, legendrecoeffs)
  void FORT_LEGENDRECOEFFS
  (
   Real* coeffs,
   const int* l,
   const int* m
   );

#define FORT_LEGENDREEVAL FORTRAN_NAME(LEGENDREEVAL, legendreeval)
  Real FORT_LEGENDREEVAL
  (
   const Real* x,
   Real* coeffs,
   const int* l,
   const int* m
   );

#define FORT_LEGENDREDERIV FORTRAN_NAME(LEGENDREDERIV, legendrederiv)
  Real FORT_LEGENDREDERIV
  (
   const Real* x,
   Real* coeffs,
   const int* l,
   const int* m
   );
}

#endif
