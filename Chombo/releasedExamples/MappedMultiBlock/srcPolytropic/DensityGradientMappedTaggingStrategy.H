#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DENSITYGRADIENTMAPPEDTAGGINGSTRATEGY_H_
#define _DENSITYGRADIENTMAPPEDTAGGINGSTRATEGY_H_

#include "AMRLevelMappedTaggingStrategy.H"

#include "NamespaceHeader.H"

//! \class DensityGradientMappedTaggingStrategy
//! This class implements a tagging strategy for the advection of a
//! DensityGradient profile.
class DensityGradientMappedTaggingStrategy: public AMRLevelMappedTaggingStrategy
{
  public:

  //! Constructor.
  DensityGradientMappedTaggingStrategy(Real a_refinementThreshold,
                                       bool a_refinementIsScaled,
                                       bool a_tagPressure,
                                       bool a_tagVorticity);

  /// Destructor.
  ~DensityGradientMappedTaggingStrategy();

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

  /// Compute max vorticity of the state, and store their magnitudes in a_vecMag.
  Real maxVorticity(LevelData<FArrayBox>& a_vecMag,
                    LevelGridMetrics* a_gridMetrics,
                    const MOLPhysics& a_physics,
                    const MappedLevelData& a_data,
                    Real a_dx) const;

  private:

  Real m_refineThresh;
  bool m_refinementIsScaled, m_tagPressure, m_tagVorticity;

};

#include "NamespaceFooter.H"

#endif
