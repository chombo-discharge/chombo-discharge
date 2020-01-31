#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDTAGGINGSTRATEGY_H_
#define _AMRLEVELMAPPEDTAGGINGSTRATEGY_H_

#include "IntVectSet.H"
#include "MOLPhysics.H"
#include "LevelGridMetrics.H"
#include "MappedLevelData.H"

#include "NamespaceHeader.H"

//! \class AMRLevelMappedTaggingStrategy
//! This class defines an interface for tagging cells for refinement
//! in high-order mapped AMR algorithms.
class AMRLevelMappedTaggingStrategy
{
  public:

  //! Base class strategy constructor. Must be called by all subclasses.
  AMRLevelMappedTaggingStrategy();

  /// Destructor.
  virtual ~AMRLevelMappedTaggingStrategy();

  //! Override this method to return the set of cells to be refined
  //! on the computational domain.
  virtual IntVectSet cellsToRefine(int a_level,
                                   const ProblemDomain& a_problem_domain,
                                   const DisjointBoxLayout& a_grids,
                                   LevelGridMetrics* a_gridMetrics,
                                   const MOLPhysics& a_physics,
                                   const MappedLevelData& a_data,
                                   Real a_dx,
                                   Real a_time) const = 0;

  //! Override this method to create a new AMRLevelMappedStabilityStrategy instance.
  virtual AMRLevelMappedTaggingStrategy* clone() const = 0;

  private:

  // Copy and assignment not allowed
  AMRLevelMappedTaggingStrategy(const AMRLevelMappedTaggingStrategy&);
  AMRLevelMappedTaggingStrategy& operator=(const AMRLevelMappedTaggingStrategy&);
};

#include "NamespaceFooter.H"

#endif
