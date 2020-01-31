#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if 0

// DTGraves, March, 2000

// an interface for EB  data to Vince Beckner's ArrayView
// include this file in application code

#ifndef _EBARRAYVIEW_H_
#define _EBARRAYVIEW_H_

#include "EBArrayViewClient.H"
#include "NamespaceHeader.H"

//the name says it all.
//this function is just to make sure that
//the functions are in your namespace
//in the debugger.
void EBArrayViewNeverShouldBeCalled()
{
  EBCellFABView(NULL);
  EBFaceFABView(NULL);
  IVFABView(NULL);
  IFFABView(NULL);
  MultiEBCellFABView(NULL);
  LDEBCellFABView(NULL, NULL);
}
#include "NamespaceFooter.H"
#endif

#endif // 0
