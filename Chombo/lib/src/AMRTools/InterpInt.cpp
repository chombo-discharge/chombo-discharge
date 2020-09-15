#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "InterpInt.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
InterpInt::
InterpInt(const DisjointBoxLayout& a_gridsFine,
          const DisjointBoxLayout& a_gridsCoar,
          int a_ref_ratio,
          int a_num_comps)
{
  define(a_gridsFine, a_gridsCoar, a_ref_ratio, a_num_comps);
}


void 
InterpInt::
define(const DisjointBoxLayout& a_gridsFine,
       const DisjointBoxLayout& a_gridsCoar,
       int a_ref_ratio,
       int a_num_comps)
{
  m_isDefined = true;
  m_refRat = a_ref_ratio;
  m_nComp  = a_num_comps;
  m_dblFine = a_gridsFine;
  m_dblCoar = a_gridsCoar;
  coarsen(m_dblCoFi, a_gridsFine, m_refRat);
  m_bufCoFi.define(m_dblCoFi, a_num_comps, IntVect::Zero);
}


///
void
InterpInt::
fillFineFromCoarPWConst(      LevelData<BaseFab<int> > & a_fine,
                        const LevelData<BaseFab<int> > & a_coar,
                        int a_srcComp,
                        int a_dstComp,
                        int a_nComp)
{
  CH_assert(m_isDefined);
  CH_assert(a_nComp <= m_nComp);
  Interval srcComps(a_srcComp, a_srcComp + a_nComp -1);
  Interval dstComps(0, a_nComp-1);
  a_coar.copyTo(srcComps, m_bufCoFi, dstComps);
  
  for(DataIterator dit = m_dblFine.dataIterator(); dit.ok(); ++dit)
  {
    const Box& fineBox = m_dblFine[dit()];
    for(BoxIterator bit(fineBox); bit.ok(); ++bit)
    {
      const IntVect& fineIV = bit();
      IntVect coarIV= coarsen(fineIV, m_refRat);
      for(int icomp = 0; icomp < a_nComp; icomp++)
      {
        int idst = a_dstComp + icomp;
        int isrc = a_srcComp + icomp;
        
        a_fine[dit()](fineIV, idst) = m_bufCoFi[dit()](coarIV, isrc);
      }
    }
  }
}


#include "NamespaceFooter.H"
