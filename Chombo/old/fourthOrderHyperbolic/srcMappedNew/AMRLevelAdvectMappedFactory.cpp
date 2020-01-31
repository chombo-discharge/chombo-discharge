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
#include "AMRLevelAdvectMapped.H"
#include "AMRLevelAdvectMappedFactory.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
AMRLevelAdvectMappedFactory::AMRLevelAdvectMappedFactory()
{
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////
// Virtual constructor
AMRLevel* AMRLevelAdvectMappedFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelAdvectMapped
  AMRLevelAdvectMapped* amrAdvectPtr = new AMRLevelAdvectMapped();

  // Set up new object
  transferSettings(amrAdvectPtr);
  // Return it
  return (static_cast <AMRLevel*> (amrAdvectPtr));
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelAdvectMappedFactory::~AMRLevelAdvectMappedFactory()
{
}

#include "NamespaceFooter.H"
