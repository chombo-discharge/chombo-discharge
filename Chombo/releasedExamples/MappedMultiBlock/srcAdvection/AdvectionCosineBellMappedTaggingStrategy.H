#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONCOSINEBELLMAPPEDTAGGINGSTRATEGY_H_
#define _ADVECTIONCOSINEBELLMAPPEDTAGGINGSTRATEGY_H_

#include "AMRLevelMappedTaggingStrategy.H"

#include "NamespaceHeader.H"

//! \class AdvectionCosineBellMappedTaggingStrategy
//! This class implements a tagging strategy for the advection of a
//! AdvectionCosineBell profile.
class AdvectionCosineBellMappedTaggingStrategy: public AMRLevelMappedTaggingStrategy
{
public:

  //! Constructor.
  AdvectionCosineBellMappedTaggingStrategy(
                                         Real a_ambientDensity,
                                         Real a_deltaDensity,
                                         Real a_size,
                                         Real a_advectVelocity,
                                         Real a_advectAngle,
                                         RealVect a_center,
                                         Real a_refinementThreshold);

  /// Destructor.
  ~AdvectionCosineBellMappedTaggingStrategy();

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

  Real m_ambientDensity;

  Real m_deltaDensity;

  Real m_size;

  Real m_advectVelocity;

  Real m_advectAngle;

  RealVect m_center;

  Real m_refineThresh;
};

#include "NamespaceFooter.H"

#endif
