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
void AMRLevelAdvectMapped::setDataMapped(LevelData<FArrayBox>& a_U,
                                         Real a_time,
                                         bool a_includeJ) const
{
  AMRLevelMappedCons::setDataMapped(a_U, a_time, a_includeJ);

  LevelAdvectMappedOperator& levelAdvectOperator =
    (LevelAdvectMappedOperator&) *m_levelConsOperatorPtr;
  levelAdvectOperator.setAdvVel(0.);
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelAdvectMapped::computeNewDt()
{
  LevelAdvectMappedOperator& levelAdvectOperator =
    (LevelAdvectMappedOperator&) *m_levelConsOperatorPtr;
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      // petermc, 29 Oct 2010, multiplied by m_cfl
      LevelData<FArrayBox> cellAvgU(m_grids, m_numStates, m_ghostVect);
      // cellAvgU = m_Unew / cellAvgJ
      levelAdvectOperator.cellUJToCellU(cellAvgU, m_Unew);
      newDT = m_cfl * m_dx / getMaxWaveSpeed(cellAvgU);
    }
  else
    {
      LevelData<FluxBox>& advVelFace = levelAdvectOperator.advVelFace();
      newDT = computeMappedNewDt(advVelFace, m_coordSysPtr, m_cfl);
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
