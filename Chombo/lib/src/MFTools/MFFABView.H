#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFFABVIEW_H_
#define _MFFABVIEW_H_

#include "MFAMRIO.H"
#include "NamespaceHeader.H"

int InitWriteMFFAB()
{
#ifdef CH_USE_HDF5
  writeMFFAB(NULL, 0);

  viewMFFAB(NULL, 0);

  writeMFFABname(NULL, 0, "label");

  writeMFLevel(NULL, 0);

  viewMFLevel(NULL, 0);

  writeMFLevelname(NULL, 0, "label");

  writeMFAMR(NULL, 0);

  viewMFAMR(NULL, 0);

  writeMFAMRname(NULL, 0, "label");
#endif
  return 1;
}

#include "NamespaceFooter.H"
#endif
