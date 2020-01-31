#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhysMappedIBC.H"
#include "SingleBlockCSAdaptor.H"
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysMappedIBC::PhysMappedIBC()
  :
  m_coordSysPtr(NULL),
  m_time(0.),
  m_haveCoordSys(false),
  m_haveTime(false),
  m_haveExactSoln(false)
{
  m_isDefined = false;
}

PhysMappedIBC::~PhysMappedIBC()
{
}

// Define the object
void PhysMappedIBC::define(const ProblemDomain& a_domain,
                           const Real&          a_dx)
{
  m_domain = a_domain;
  m_dx     = a_dx;

  m_isDefined = true;
}

void PhysMappedIBC::setTime(Real a_time)
{
  m_time = a_time;
  m_haveTime = true;
}

void PhysMappedIBC::setCoordSys(MultiBlockCoordSys* a_coordSysPtr)
{
  CH_assert(a_coordSysPtr != NULL);
  m_coordSysPtr = a_coordSysPtr;
  m_haveCoordSys = true;
}

bool PhysMappedIBC::haveExactSoln() const
{
  return m_haveExactSoln;
}

//-----------------------------------------------------------------------
void
PhysMappedIBC::
setCoordSys(NewCoordSys* a_coordSysPtr)
{
  MultiBlockCoordSys* coordSys = new SingleBlockCSAdaptor(a_coordSysPtr, m_domain);
  setCoordSys(coordSys);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
