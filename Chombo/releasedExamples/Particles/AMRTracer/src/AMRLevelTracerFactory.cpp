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

#include "AMRLevelTracerFactory.H"

#include "NamespaceHeader.H"

AMRLevelTracerFactory::
AMRLevelTracerFactory(AdvectionVelocityFunction   a_advFunc,
                      const Real&                   a_cfl,
                      const Real&                 a_domainLength,
                      const int&                  a_maxGridSize)
{
  m_advFunc                = a_advFunc;
  m_cfl                    = a_cfl;
  m_domainLength           = a_domainLength;
  m_maxGridSize            = a_maxGridSize;
}

// Virtual constructor
AMRLevel* AMRLevelTracerFactory::new_amrlevel() const
{
  // Create a new AMRLevelTracer
  AMRLevelTracer* amrPartPtr = new AMRLevelTracer(m_advFunc,
                                                  m_cfl,
                                                  m_domainLength,
                                                  m_maxGridSize);

  // Return it
  return (static_cast <AMRLevel*> (amrPartPtr));
}

#include "NamespaceFooter.H"
