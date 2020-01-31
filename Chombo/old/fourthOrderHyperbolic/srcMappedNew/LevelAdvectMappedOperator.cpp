#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelAdvectMappedOperator.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelAdvectMappedOperator::LevelAdvectMappedOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchAdvectMappedOperator();
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelAdvectMappedOperator::~LevelAdvectMappedOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
void LevelAdvectMappedOperator::setPatchIndex(const DataIndex&  a_ind) const
{
  // This function must be const because cellUJToCellU must be const.
  // ugh, need to do this for the function to be const
  PatchAdvectMappedOperator& patchConsOperator =
    (PatchAdvectMappedOperator&) *m_patchConsOperatorPtr;

  const Box& curBox = m_grids[a_ind];
  patchConsOperator.setCurrentBox(curBox);
  patchConsOperator.setAdvVelAvg(&(m_advVelAvg[a_ind]));
  patchConsOperator.setAdvVelCen(&(m_advVelCen[a_ind]));
}

//////////////////////////////////////////////////////////////////////////////
void
LevelAdvectMappedOperator::setAdvVel(const LevelData<FArrayBox>&  a_advVel)
{
  const IntVect& ghostVect = a_advVel.ghostVect();
  // a_advVel is cell-centered
  // m_advVelCen is face-centered
  m_advVelCen.define(m_grids, SpaceDim, ghostVect);
  computeFaceAverages(m_advVelCen, a_advVel);

  // m_advVelAvg is face-averaged
  m_advVelAvg.define(m_grids, SpaceDim, ghostVect);

  PatchAdvectMappedOperator& patchAdvectOperator =
    (PatchAdvectMappedOperator&) *m_patchConsOperatorPtr;

  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      patchAdvectOperator.deconvolveAdvVel();
    }
}

#include "NamespaceFooter.H"
