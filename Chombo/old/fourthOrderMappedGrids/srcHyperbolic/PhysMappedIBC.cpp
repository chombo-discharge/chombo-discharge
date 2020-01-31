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
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysMappedIBC::PhysMappedIBC()
{
  m_isDefined = false;

  m_coordSysPtr = NULL;
  m_time = 0.;
  m_gotTime = false;
  m_gotCoordSys = false;
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
  m_gotTime = true;
}


void PhysMappedIBC::setCoordSys(MultiBlockCoordSys* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  m_gotCoordSys = true;
}

#include "NamespaceFooter.H"
