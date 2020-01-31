#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELADVECTMAPPEDFACTORY_H_
#define _AMRLEVELADVECTMAPPEDFACTORY_H_

#include "AMRLevelMappedConsFactory.H"
#include "AMRLevelAdvectMapped.H"

#include "NamespaceHeader.H"

/// AMR Advection Equation factory
/**
 */
class AMRLevelAdvectMappedFactory : public AMRLevelMappedConsFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelAdvectMappedFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelAdvectMappedFactory();
};

#include "NamespaceFooter.H"

#endif
