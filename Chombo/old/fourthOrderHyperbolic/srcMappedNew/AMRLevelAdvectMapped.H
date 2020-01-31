#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELADVECTMAPPED_H_
#define _AMRLEVELADVECTMAPPED_H_

#include "AMRLevelMappedCons.H"

#include "NamespaceHeader.H"

/// AMR Level for fourth-order solution to conservation law
/**
 */
class AMRLevelAdvectMapped : public AMRLevelMappedCons
{
  friend class AMRLevelAdvectMappedFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelAdvectMapped();

  /// Destructor
  /**
   */
  virtual ~AMRLevelAdvectMapped();

  /// Initialize data
  virtual void initialData();

protected:

  /// Compute dt using m_Unew
  virtual Real computeNewDt();

  /// Index within primitive variables for tagging cells
  virtual int indexForTagging();

  // cell-centered advection velocity
  LevelData<FArrayBox> m_advVel;
};

#include "NamespaceFooter.H"
#endif
