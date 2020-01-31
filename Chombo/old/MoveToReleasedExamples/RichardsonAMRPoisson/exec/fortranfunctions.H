#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FORTRANFUNCTIONS_H_
#define _FORTRANFUNCTIONS_H_

#include "FORT_PROTO.H"
#include "REAL.H"

extern "C"
{

// Prototype for Fortran procedure setrhs ...
//
#define FORT_SETRHS FORTRAN_NAME( SETRHS ,setrhs )
void
FORT_SETRHS(
      CHFp_FRA(rhs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(dx) );


#define FORT_SETEXACT FORTRAN_NAME( SETEXACT, setexact )
void
FORT_SETEXACT(
      CHFp_FRA(rhs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(dx) );
}
#endif
