#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFDEBUGOUT_H_
#define _MFDEBUGOUT_H_

#include "EBDebugOut.H"
#include "NamespaceHeader.H"

class MFCellFAB;

extern "C"
{
  ///
  void
  dumpLDMFCF(const LevelData<MFCellFAB>*  ldptr);

  ///
  void
  getMaxMFLevel(const LevelData<MFCellFAB>*  ldptr);

  ///
  void
  getMaxMFFAB(const MFCellFAB*  cfptr);


  ///
  /**
   */
  void dumpMFLDDBL(const LevelData<MFCellFAB>* a_stenPtr);

  ///
  /**
     Dump inputs to standard out.  For use inside debuggers.
   */
  void dumpLevDBIVF(const LevelData< BaseIVFAB<Real> >* a_ldPtr);

}

#include "NamespaceFooter.H"
#endif
