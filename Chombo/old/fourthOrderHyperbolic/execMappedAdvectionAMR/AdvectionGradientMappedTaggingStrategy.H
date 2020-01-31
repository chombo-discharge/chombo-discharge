#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONGRADIENTMAPPEDTAGGINGSTRATEGY_H_
#define _ADVECTIONGRADIENTMAPPEDTAGGINGSTRATEGY_H_

#include "AMRLevelMappedTaggingStrategy.H"

#include "NamespaceHeader.H"

//! \class AdvectionGradientMappedTaggingStrategy
//! This class implements a tagging strategy for the advection of a
//! AdvectionGradient profile.
class AdvectionGradientMappedTaggingStrategy: public AMRLevelMappedTaggingStrategy
{
  public:

  //! Constructor.
  AdvectionGradientMappedTaggingStrategy(Real a_refinementThreshold,
                                         bool a_refinementIsScaled,
                                         bool a_relativeGradient);

  /// Destructor.
  ~AdvectionGradientMappedTaggingStrategy();

  // Overriden methods

  IntVectSet cellsToRefine(int a_level,
                           const ProblemDomain& a_problem_domain,
                           const DisjointBoxLayout& a_grids,
                           LevelGridMetrics* a_gridMetrics,
                           const MOLPhysics& a_physics,
                           const MappedLevelData& a_data,
                           Real a_dx,
                           Real a_time) const;

  AMRLevelMappedTaggingStrategy* clone() const;

  private:

  Real m_refineThresh;
  bool m_refinementIsScaled;
  bool m_relativeGradient;
};

#include "NamespaceFooter.H"

#endif
