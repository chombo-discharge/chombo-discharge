#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ANALYTICLONLATADVECTIONMAPPEDSTABILITYSTRATEGY_H_
#define _ANALYTICLONLATADVECTIONMAPPEDSTABILITYSTRATEGY_H_

#include "AMRLevelMappedStabilityStrategy.H"
#include "RefCountedPtr.H"
#include "VectorFunction.H"
#include "ScalarFunction.H"

#include "NamespaceHeader.H"

//! \class AnalyticLonLatAdvectionMappedStabilityStrategy
//! This class implements a stability strategy for the Euler equations.
//! in high-order mapped AMR algorithms.
class AnalyticLonLatAdvectionMappedStabilityStrategy: public AMRLevelMappedStabilityStrategy
{
  public:

  //! Constructor.
  //! \param a_stabilityFactor The factor representing the broadest height of the region
  //!                          of stability for the time integrator of interest.
  AnalyticLonLatAdvectionMappedStabilityStrategy(Real a_stabilityFactor,
                                                 RefCountedPtr<VectorFunction> a_velocity);

  //! Stream function versioin of above constructor.
  //! \param a_stabilityFactor The factor representing the broadest height of the region
  //!                          of stability for the time integrator of interest.
  AnalyticLonLatAdvectionMappedStabilityStrategy(Real a_stabilityFactor,
                                                 RefCountedPtr<ScalarFunction> a_streamFunc);

  /// Destructor.
  ~AnalyticLonLatAdvectionMappedStabilityStrategy();

  // Overriden methods.

  Real maxWaveSpeed(const ProblemDomain& a_problem_domain,
                    const DisjointBoxLayout& a_grids,
                    const LevelGridMetrics* const a_gridMetrics,
                    const MOLPhysics& a_physics,
                    const MappedLevelData& a_data,
                    Real a_dx,
                    Real a_time) const;

  AMRLevelMappedStabilityStrategy* clone() const;

  private:

  RefCountedPtr<VectorFunction> m_velocity;
  RefCountedPtr<ScalarFunction> m_streamFunc;

};

#include "NamespaceFooter.H"

#endif
