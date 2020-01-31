#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeArrayView.H

#ifndef _NODEARRAYVIEW_H_
#define _NODEARRAYVIEW_H_

#include "ArrayView.H"
#include "NodeDatasetClient.H"
#include "NamespaceHeader.H"


#ifndef NDEBUG
// Force these debug functions to be linked into executables
void NodeNeverCalled()
{
  NeverCalled();

  // For use with NodeFArrayBox
  ArrayViewNode(NULL);

  ArrayViewBool(NULL);

  MultiArrayViewNodeFab(NULL);

  MultiArrayViewLevelNodeFab(NULL);

  MultiArrayViewInt(NULL);

}
#endif

#include "NamespaceFooter.H"
#endif
