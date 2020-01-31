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
void
LevelAdvectMappedOperator::define(
                                  const DisjointBoxLayout&  a_thisLayout,
                                  const DisjointBoxLayout&  a_coarserLayout,
                                  const ProblemDomain&      a_domain,
                                  const int&                a_refineCoarse,
                                  const Real&               a_dx,
                                  const MOLPhysics* const a_molPhysics,
                                  const int&                a_numFields,
                                  const bool&               a_hasCoarser,
                                  const bool&               a_hasFiner)
{
  LevelConsOperator::define(a_thisLayout, a_coarserLayout,
                            a_domain, a_refineCoarse, a_dx, a_molPhysics,
                            a_numFields, a_hasCoarser, a_hasFiner);
  m_physIBCPtr = (PhysAdvectMappedIBC*) a_molPhysics->getPhysIBC();
  defineAdvVel(); // after having set m_grids and m_numGhost
}


//////////////////////////////////////////////////////////////////////////////
void
LevelAdvectMappedOperator::defineAdvVel()
{
  IntVect ghostVect = m_numGhost * IntVect::Unit;
  // cell-centered LevelData<FArrayBox>
  m_advVel.define(m_grids, SpaceDim, ghostVect);
  // face-centered LevelData<FluxBox>
  m_advVelFace.define(m_grids, SpaceDim, ghostVect);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelAdvectMappedOperator::evalRHS(
                                   LevelData<FArrayBox>&       a_LofU,
                                   LevelData<FArrayBox>&       a_U,
                                   LevelFluxRegister&          a_finerFluxRegister,
                                   LevelFluxRegister&          a_coarserFluxRegister,
                                   const LevelData<FArrayBox>& a_UcoarseOld,
                                   const Real&                 a_timeCoarseOld,
                                   const LevelData<FArrayBox>& a_UcoarseNew,
                                   const Real&                 a_timeCoarseNew,
                                   Real                        a_time,
                                   Real                        a_weight)
{
  setAdvVel(a_time);
  // Now that we have advection velocities, do what we did before.
  LevelMappedConsOperator::evalRHS(a_LofU, a_U,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   a_UcoarseOld, a_timeCoarseOld,
                                   a_UcoarseNew, a_timeCoarseNew,
                                   a_time, a_weight);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelAdvectMappedOperator::setAdvVel(Real a_time)
{
  // These lines added 29 March 2011, so that restart will work properly.
  m_physIBCPtr->setCoordSys(m_coordSysPtr);
  m_physIBCPtr->setTime(a_time);
  // Set cell-centered advection velocity.
  m_physIBCPtr->advVel(m_advVel, a_time);
  // Use m_advVel to set face-centered m_advVelFace
  IntVect ghostVect = m_numGhost * IntVect::Unit;
  // m_advVel is cell-centered, and we now have it.
  // m_advVelFace is face-centered; calculate it from cell-centered m_advVel.
  // m_advVelFace.define(m_grids, SpaceDim, ghostVect);
  computeFaceCenters(m_advVelFace, m_advVel);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelAdvectMappedOperator::setPatchIndex(const DataIndex&  a_ind) const
{
  // This function must be const because cellUJToCellU must be const.
  // ugh, need to do this for the function to be const
  PatchAdvectMappedOperator& patchConsOperator =
    (PatchAdvectMappedOperator&) *m_patchConsOperatorPtr;

  const Box& curBox = m_grids[a_ind];
  patchConsOperator.setCurrentBox(curBox);
  patchConsOperator.setAdvVelFace(&(m_advVelFace[a_ind]));
}

#include "NamespaceFooter.H"
