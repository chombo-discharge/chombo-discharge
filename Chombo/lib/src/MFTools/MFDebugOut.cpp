#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFDebugOut.H"
#include "MFCellFAB.H"
#include "EBDebugOut.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "FaceIndex.H"
#include "Vector.H"
#include "VolIndex.H"
#include "LoHiSide.H"
#include "VoFIterator.H"
#include "DebugOut.H"
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include "NamespaceHeader.H"

void
getMaxMFFAB(const MFCellFAB*  cfptr)
{
  int numPhase = cfptr->numPhases();
  for (int iphase=0; iphase<numPhase; iphase++)
    {
      const EBCellFAB& fab = cfptr->getPhase(iphase);
      for (int icomp = 0; icomp < fab.nComp(); icomp++)
        {
          Real maxval = -1.0e16;
          Real minval =  1.0e16;
          VolIndex vofmax, vofmin;
          pout() << "phase = " << iphase << "; c = " << icomp << ", ";

          const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
          const Box& grid = fab.box();
          IntVectSet ivs(grid);
          for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
            {
              Real datval = fab(vofit(), icomp);
              if (datval > maxval)
                {
                  maxval = datval;
                  vofmax = vofit();
                }
              if (datval < minval)
                {
                  minval = datval;
                  vofmin = vofit();
                }
            }
          pout() << "max=" <<  maxval << " at " << vofmax << ", ";
          pout() << "min=" <<  minval << " at " << vofmin << endl;
        }
    }
}

void
getMaxMFLevel(const LevelData<MFCellFAB>*  ldptr)
{
  int numPhase = 2;
  const LevelData<MFCellFAB>& ld = *ldptr;
  for (int iphase=0; iphase<numPhase; iphase++)
    {
      for (int icomp = 0; icomp < ld.nComp(); icomp++)
        {
          Real maxval = -1.0e16;
          Real minval =  1.0e16;
          VolIndex vofmax, vofmin;
          pout() << "phase = " << iphase << "; c = " << icomp << ", ";
          const DisjointBoxLayout& dbl = ld.disjointBoxLayout();
          for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
            {
              const EBCellFAB& fab = ld[dit()].getPhase(iphase);
              const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
              const Box& grid = dbl.get(dit());
              IntVectSet ivs(grid);
              for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
                {
                  Real datval = fab(vofit(), icomp);
                  if (datval > maxval)
                    {
                      maxval = datval;
                      vofmax = vofit();
                    }
                  if (datval < minval)
                    {
                      minval = datval;
                      vofmin = vofit();
                    }
                }
            }
          pout() << "max=" <<  maxval << " at " << vofmax << ", ";
          pout() << "min=" <<  minval << " at " << vofmin << endl;
        }
    }
}
void
dumpLDMFCF(const LevelData<MFCellFAB>*  ldptr)
{
  int numPhase = 2;
  const LevelData<MFCellFAB>& ld = *ldptr;
  for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      for (int iphase=0; iphase<numPhase; iphase++)
        {
          const FArrayBox& singFAB = (FArrayBox&)ld[dit()].getPhase(iphase).getSingleValuedFAB();
          const MiniIVFAB<Real>& multFAB = ld[dit()].getPhase(iphase).getMultiValuedFAB();
          dumpFAB(&singFAB);
          dumpIVFAB(&multFAB);
        }
    }
}

void
dumpMFLDDBL(const LevelData<MFCellFAB>*  memLDF_Ptr)
{
  const DisjointBoxLayout& dbl = memLDF_Ptr->disjointBoxLayout();
  dumpDBL(&dbl);
}

void dumpLevDBIVF(const LevelData< BaseIVFAB<Real> >* a_ldptr)
{
  const LevelData< BaseIVFAB<Real> >& ld = *a_ldptr;
  for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      dumpIVFAB(&ld[dit()]);
    }
}
#include "NamespaceFooter.H"
