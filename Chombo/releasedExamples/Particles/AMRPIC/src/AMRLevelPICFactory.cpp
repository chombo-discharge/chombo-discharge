#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"

#include "AMRLevelPICFactory.H"

#include "NamespaceHeader.H"

AMRLevelPICFactory::
AMRLevelPICFactory()
{
  m_ibc = NULL;
}

AMRLevelPICFactory::
~AMRLevelPICFactory()
{
  if (m_ibc != NULL)
  {
    delete m_ibc;
  }
}

AMRLevelPICFactory::
AMRLevelPICFactory(const ZeldovichIBC*         a_ibc,
		   const Real&                 a_cfl,
		   const Real&                 a_domainLength,
		   const int&                  a_maxGridSize,
		   const InterpType&           a_poissonInterpType,
		   const Real&                 a_refineThresh)
{

  m_ibc = a_ibc->new_zeldovichIBC();

  m_cfl                    = a_cfl;
  m_domainLength           = a_domainLength;
  m_maxGridSize            = a_maxGridSize;
  m_poissonInterpType      = a_poissonInterpType;
  m_refineThresh           = a_refineThresh;
}

// Virtual constructor
AMRLevel* AMRLevelPICFactory::new_amrlevel() const
{
  // Create a new AMRLevelPIC
  AMRLevelPIC* amrPartPtr = new AMRLevelPIC(m_ibc,
					    m_cfl,
					    m_domainLength,
					    m_maxGridSize,
					    m_poissonInterpType,
					    m_refineThresh);

  // Return it
  return (static_cast <AMRLevel*> (amrPartPtr));
}

#include "NamespaceFooter.H"
