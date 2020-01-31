#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HEATFORT_H_
#define _HEATFORT_H_

#include "REAL.H"
#include "FORT_PROTO.H"

extern "C"
{
 void testfort_(const int* const argint, const Real* const argreal);
 void testfortbigf_(const int* const argint, const Real* const argreal);
};

#endif
