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
#include "AMRLevelShallowWaterMapped.H"
#include "AMRLevelShallowWaterMappedFactory.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
AMRLevelShallowWaterMappedFactory::AMRLevelShallowWaterMappedFactory()
{
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////
// Virtual constructor
AMRLevel* AMRLevelShallowWaterMappedFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelShallowWaterMapped
  AMRLevelShallowWaterMapped* amrShallowWaterPtr =
      new AMRLevelShallowWaterMapped();

  // Set up new object
  transferSettings(amrShallowWaterPtr);
  // Return it
  return (static_cast <AMRLevel*> (amrShallowWaterPtr));
}

//////////////////////////////////////////////////////////////////////////////

AMRLevelShallowWaterMappedFactory::~AMRLevelShallowWaterMappedFactory()
{
}

#include "NamespaceFooter.H"
