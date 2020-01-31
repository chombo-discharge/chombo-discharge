#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELSHALLOWWATERMAPPEDFACTORY_H_
#define _AMRLEVELSHALLOWWATERMAPPEDFACTORY_H_

#include "AMRLevelMappedConsFactory.H"
#include "AMRLevelShallowWaterMapped.H"

#include "NamespaceHeader.H"

/// AMR Shallow Water Equation factory
/**
 */
class AMRLevelShallowWaterMappedFactory : public AMRLevelMappedConsFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelShallowWaterMappedFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelShallowWaterMappedFactory();
};

#include "NamespaceFooter.H"

#endif
