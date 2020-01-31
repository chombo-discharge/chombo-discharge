#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectConstantIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

// Null constructor
AdvectConstantIBC::AdvectConstantIBC()
{
  m_haveAdvVel = false;
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
AdvectConstantIBC::AdvectConstantIBC(
                                     Real a_fieldVal,
                                     Real a_advectVelocity,
                                     Real a_advectAngle)
{
  m_fieldVal = a_fieldVal;
  m_advectVelocity = a_advectVelocity;
  m_advectAngle = a_advectAngle;
  m_haveAdvVel = false;
}

AdvectConstantIBC::~AdvectConstantIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* AdvectConstantIBC::new_physIBC()
{
  AdvectConstantIBC* retval = new AdvectConstantIBC();
  retval->m_fieldVal = m_fieldVal;
  retval->m_advectVelocity = m_advectVelocity;
  retval->m_advectAngle = m_advectAngle;

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void AdvectConstantIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  // int nComp = a_U.nComp();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Storage for current grid
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(m_fieldVal);
    }
}

//-----------------------------------------------------------------------

// Return the advection velocity at each point
void AdvectConstantIBC::advVel(FArrayBox& a_advVelFab,
                               Real a_time)
{
  CH_TIME("AdvectConstantIBC::advVel");

  a_advVelFab.setVal(m_fieldVal);
  // m_haveAdvVel = true;

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

#include "NamespaceFooter.H"
