#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VOLTESTPARAMETERS_H_
#define _VOLTESTPARAMETERS_H_

#include "REAL.H"

namespace VT
{

enum ICType
{
  ICTypeConstant,
  ICTypeLinear
};

enum RefType
{
  RefTypeCenter,
  RefTypeRotating
};

extern int g_runStatus;
extern Real g_conservedJU0;
extern Real g_conservedJU1;

}  // End of namespace VT

#endif
