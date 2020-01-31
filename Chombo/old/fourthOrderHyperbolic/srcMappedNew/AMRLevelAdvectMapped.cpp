#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include <string>
#include "parstream.H"

#include "AMRLevelAdvectMapped.H"
#include "LevelAdvectMappedOperator.H"

#include "SPMD.H"
#include "computeMappedNewDt.H"
// #include "DebugInclude.H"
// for kludge
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelAdvectMapped::AMRLevelAdvectMapped()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvectMapped default constructor" << endl;
  }
  m_levelConsOperatorPtr = new LevelAdvectMappedOperator();
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
AMRLevelAdvectMapped::~AMRLevelAdvectMapped()
{
}

//////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelAdvectMapped::initialData()
{
  AMRLevelMappedCons::initialData();

  // need to define m_advVel
  m_advVel.define(m_grids, SpaceDim, m_ghostVect);
  PhysAdvectMappedIBC* physIBCPtr =
    (PhysAdvectMappedIBC*) m_molPhysics->getPhysIBC();
  physIBCPtr->advVel(m_advVel);

  LevelAdvectMappedOperator& levelAdvectOperator =
    (LevelAdvectMappedOperator&) *m_levelConsOperatorPtr;

  levelAdvectOperator.setAdvVel(m_advVel);
  // need to define m_advVelAvg and m_advVelCen in LevelAdvectMappedOperator
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelAdvectMapped::computeNewDt()
{
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      // petermc, 29 Oct 2010, multiplied by m_cfl
      newDT = m_cfl * m_dx / getMaxWaveSpeed(m_Unew);
    }
  else
    {
      LevelAdvectMappedOperator& levelAdvectOperator =
        (LevelAdvectMappedOperator&) *m_levelConsOperatorPtr;
      LevelData<FluxBox>& advVelAvg = levelAdvectOperator.advVelAvg();
      newDT = computeMappedNewDt(advVelAvg, m_coordSysPtr, m_cfl);
    }
  m_dtNew = newDT;
  return newDT;
}

//////////////////////////////////////////////////////////////////////////////

int AMRLevelAdvectMapped::indexForTagging()
{
  // tag on the advected quantity
  return 0;
}

#include "NamespaceFooter.H"
