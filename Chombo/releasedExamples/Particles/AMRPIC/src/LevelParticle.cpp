#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelParticle.H"

#include "NamespaceHeader.H"

// Constructor - set up some defaults
LevelParticle::LevelParticle()
{
  m_defined = false;
  m_dx = 0.0;
  m_patchParticle = NULL;
}

// Destructor - free up storage
LevelParticle::~LevelParticle()
{
  if (m_patchParticle != NULL)
    {
      delete m_patchParticle;
    }
}

// Define the object so that time stepping can begin
void LevelParticle::define(const ProblemDomain&     a_domain,
                           const Real&              a_dx,
                           const PatchParticle* const  a_patchParticleFactory)
{
  // Sanity checks
  CH_assert(a_dx>0);

  // Cache data
  m_dx     = a_dx;
  m_domain = a_domain;

  m_domainLength = RealVect(D_DECL(a_dx*a_domain.size(0),
                                   a_dx*a_domain.size(1),
                                   a_dx*a_domain.size(2)));

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchParticle != NULL)
    {
      delete m_patchParticle;
    }

  m_patchParticle = a_patchParticleFactory->new_patchParticle();

  // Everything is defined
  m_defined = true;
}

#include "NamespaceFooter.H"
