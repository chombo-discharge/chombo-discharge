#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFDEBUGDUMP_H_
#define _MFDEBUGDUMP_H_

#include "MFAMRIO.H"
#include "NamespaceHeader.H"

void MFDontYouDareCallMe()
{
  viewMFFAB(NULL,0);
  browseMFFAB(NULL,0);
  viewMFLevel(NULL,0);
  browseMFLevel(NULL,0);
  viewMFAMR(NULL,0);
  browseMFAMR(NULL,0);
}
#include "NamespaceFooter.H"
#endif
