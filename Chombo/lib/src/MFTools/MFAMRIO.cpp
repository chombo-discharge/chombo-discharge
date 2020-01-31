#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Thurs, Aug 16, 2001

#include <fstream>
#include <string>
using std::fstream;
using std::string;
#include <cstdio>
#include <cmath>
#include "CH_HDF5.H"
#include "EBAMRIO.H"
#include "MFAMRIO.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "PolyGeom.H"
#include "MFAliasFactory.H"
#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

void
writeMFFAB(const MFCellFAB* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  const EBCellFAB& singlePhase = a_dataPtr->getPhase(a_phase);

  writeEBFAB(&singlePhase);
}


void
viewMFFAB(const MFCellFAB* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  const EBCellFAB& singlePhase = a_dataPtr->getPhase(a_phase);

  viewEBFAB(&singlePhase);

}


void
browseMFFAB(const MFCellFAB* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  const EBCellFAB& singlePhase = a_dataPtr->getPhase(a_phase);

  browseEBFAB(&singlePhase);


}

void
writeMFFABname(const MFCellFAB* a_dataPtr,
               int a_phase,
               const char*      a_filename)
{
  if (a_dataPtr == NULL) return;

  const EBCellFAB& singlePhase = a_dataPtr->getPhase(a_phase);

  writeEBFABname(&singlePhase, a_filename);


}

void
writeMFLevel(const LevelData<MFCellFAB>* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  LevelData<EBCellFAB> singlePhase;

  aliasMF(singlePhase, a_phase, *a_dataPtr);

  writeEBLevel(&singlePhase);

}

void
viewMFLevel(const LevelData<MFCellFAB>* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  LevelData<EBCellFAB> singlePhase;
  aliasMF(singlePhase, a_phase, *a_dataPtr);
  viewEBLevel(&singlePhase);
}

void
browseMFLevel(const LevelData<MFCellFAB>* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  LevelData<EBCellFAB> singlePhase;
  aliasMF(singlePhase, a_phase, *a_dataPtr);
  browseEBLevel(&singlePhase);

}

void
writeMFLevelname(const LevelData<MFCellFAB>* a_dataPtr,
                 int a_phase,
                 const char*                 a_filename)
{
  if (a_dataPtr == NULL) return;

  LevelData<EBCellFAB> singlePhase;
  aliasMF(singlePhase, a_phase, *a_dataPtr);
  writeEBLevelname(&singlePhase, a_filename);


}

void
writeMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr, int a_phase)
{
  if (a_dataPtr == NULL) return;

  Vector<LevelData<EBCellFAB>* > singlePhaseVect(a_dataPtr->size(), NULL);
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      singlePhaseVect[lev] = new LevelData<EBCellFAB>;
    }

  aliasMF(singlePhaseVect, a_phase, *a_dataPtr);

  writeEBAMR(&singlePhaseVect);

  // now clean up temp storage
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      delete singlePhaseVect[lev];
    }
}

void
viewMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr,
          int a_phase)
{
  if (a_dataPtr == NULL) return;

  Vector<LevelData<EBCellFAB>* > singlePhaseVect(a_dataPtr->size(), NULL);
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      singlePhaseVect[lev] = new LevelData<EBCellFAB>;
    }

  aliasMF(singlePhaseVect, a_phase, *a_dataPtr);

  viewEBAMR(&singlePhaseVect);

  // now clean up temp storage
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      delete singlePhaseVect[lev];
    }

}

void
browseMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr,
            int a_phase)
{
  if (a_dataPtr == NULL) return;


  Vector<LevelData<EBCellFAB>* > singlePhaseVect(a_dataPtr->size(), NULL);
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      singlePhaseVect[lev] = new LevelData<EBCellFAB>;
    }

  aliasMF(singlePhaseVect, a_phase, *a_dataPtr);

  browseEBAMR(&singlePhaseVect);

  // now clean up temp storage
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      delete singlePhaseVect[lev];
    }

}

void
writeMFAMRname(const Vector<LevelData<MFCellFAB>* >* a_dataPtr,
               int a_phase,
               const char*                           a_filename)
{
  if (a_dataPtr == NULL) return;

  Vector<LevelData<EBCellFAB>* > singlePhaseVect(a_dataPtr->size(), NULL);
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      singlePhaseVect[lev] = new LevelData<EBCellFAB>;
    }

  aliasMF(singlePhaseVect, a_phase, *a_dataPtr);

  writeEBAMRname(&singlePhaseVect,a_filename);

  // now clean up temp storage
  for (int lev=0; lev<singlePhaseVect.size(); lev++)
    {
      delete singlePhaseVect[lev];
    }


}


#endif   // CH_USE_HDF5
#include "NamespaceFooter.H"
