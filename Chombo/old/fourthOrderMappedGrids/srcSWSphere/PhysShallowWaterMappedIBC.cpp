#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhysShallowWaterMappedIBC.H"
#include "ShallowWaterPhysicsF_F.H"
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysShallowWaterMappedIBC::PhysShallowWaterMappedIBC() :
  PhysMappedIBC(),
  m_isFortranCommonSet(false)
{
}

PhysShallowWaterMappedIBC::~PhysShallowWaterMappedIBC()
{
}

void PhysShallowWaterMappedIBC::setFortranCommon(const Real& a_gravity,
                                            const Real& a_omega,
                                            const Real& a_alpha)
{
  CH_assert(m_isFortranCommonSet == false);

  m_gravity = a_gravity;
  m_omega = a_omega;
  m_alpha = a_alpha;

  FORT_SHALLOWWATERSETF(CHF_CONST_REAL(a_gravity),
                        CHF_CONST_REAL(a_omega),
                        CHF_CONST_REAL(a_alpha));

  m_isFortranCommonSet = true;
}

void PhysShallowWaterMappedIBC::topoHeight(LevelData<FArrayBox>& a_topoHeight)
{
}

#include "NamespaceFooter.H"
