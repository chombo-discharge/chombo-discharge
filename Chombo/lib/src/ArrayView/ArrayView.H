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

// ArrayView.H

#ifndef _ARRAYVIEW_H_
#define _ARRAYVIEW_H_

#include "DatasetClient.H"
#include "NamespaceHeader.H"


#ifndef NDEBUG
// Force these debug functions to be linked into executables
void NeverCalled()
{
  ArrayView(NULL);

  ArrayViewInt(NULL);

  MultiArrayView(NULL);

  ArrayViewFabFormatLabel(NULL, "%12.10f", "Label");

  MultiArrayViewRealBaseFab(NULL);

  MultiArrayViewFab(NULL);

  ArrayViewLDF(NULL);

  ArrayViewLDFFormatLabel(NULL, "%12.10f", "Label");

  MultiArrayViewFormatLabel(NULL, "%12.10f", "Label");

}
#endif

#include "NamespaceFooter.H"
#endif

#endif // 0
