#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// an interface for I[NE]Fab data to Vince Beckner's ArrayView
//

#ifndef _EBARRAYVIEWCLIENT_H_
#define _EBARRAYVIEWCLIENT_H_

#include "REAL.H"
#include "LayoutData.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "BaseIFFAB.H"
#include "BaseIVFAB.H"
#include "EBArrayViewClient.H"
#include "NamespaceHeader.H"
#ifdef CH_MULTIDIM
#include "EBArrayView_ExternC_Mangler.H"
#endif

extern "C"
{
  ///
  /**
     Arrayview interface for irregular data on faces.
     Can be called from debugger.
   */
  bool IFFABView(const BaseIFFAB<Real>* const fab);

  ///
  /**
     Arrayview interface for irregular data on VoFs.
     Can be called from debugger.
   */
  bool IVFABView(const BaseIVFAB<Real>* const fab);

  ///
  /**
     Arrayview interface for regular data on VoFs.
     Can be called from debugger.
   */
  bool EBCellFABView(const EBCellFAB* const fab);

  ///
  /**
     Arrayview interface for regular data on Faces.
     Can be called from debugger.
   */
  bool EBFaceFABView(const EBFaceFAB* const fab);
  bool IrregViewRealPtrArrayNVarDims(Real *data[], int nvar,
                                     const int *lodim, const int *hidim,
                                     const char *format, const char *label);

  ///
  /**
     Arrayview interface for regular data on VoFs.
     Can be called from debugger.
   */
  bool MultiEBCellFABView(const LayoutData<EBCellFAB >* const fab);

  ///
  /**
     Arrayview interface for regular data on VoFs.
     Do NOT call this one from the debugger.  Use
     only inside your program.
   */
  bool LDEBCellFABView(const LevelData<EBCellFAB >* const fab,
                       const Box* const a_domainPtr);
  void  makeFabFromEBFAB(BaseFab<Real>& fab,const EBCellFAB& ebfab);
}

#include "NamespaceFooter.H"
#endif
