#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONABSOLUTEMAPPEDTAGGINGSTRATEGY_H_
#define _ADVECTIONABSOLUTEMAPPEDTAGGINGSTRATEGY_H_

#include "AMRLevelMappedTaggingStrategy.H"

#include "NamespaceHeader.H"

//! \class AdvectionAbsoluteMappedTaggingStrategy
//! This class implements a tagging strategy for the advection of a
//! AdvectionAbsolute profile.
class AdvectionAbsoluteMappedTaggingStrategy: public AMRLevelMappedTaggingStrategy
{
public:

  //! Constructor.
  AdvectionAbsoluteMappedTaggingStrategy(Real a_refinementThreshold);

  /// Destructor.
  ~AdvectionAbsoluteMappedTaggingStrategy();

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
};

#include "NamespaceFooter.H"

#endif
