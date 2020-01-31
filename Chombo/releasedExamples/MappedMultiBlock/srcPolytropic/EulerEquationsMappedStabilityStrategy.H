#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EULEREQUATIONSMAPPEDSTABILITYSTRATEGY_H_
#define _EULEREQUATIONSMAPPEDSTABILITYSTRATEGY_H_

#include "AMRLevelMappedStabilityStrategy.H"

#include "NamespaceHeader.H"

//! \class EulerEquationsMappedStabilityStrategy
//! This class implements a stability strategy for the Euler equations.
//! in high-order mapped AMR algorithms.
class EulerEquationsMappedStabilityStrategy: public AMRLevelMappedStabilityStrategy
{
  public:

  //! Constructor.
  //! \param a_stabilityFactor The factor representing the broadest height of the region
  //!                          of stability for the time integrator of interest.
  explicit EulerEquationsMappedStabilityStrategy(Real a_stabilityFactor);

  /// Destructor.
  ~EulerEquationsMappedStabilityStrategy();

  // Overriden methods.

  Real maxWaveSpeed(const ProblemDomain& a_problem_domain,
                    const DisjointBoxLayout& a_grids,
                    const LevelGridMetrics* const a_gridMetrics,
                    const MOLPhysics& a_physics,
                    const MappedLevelData& a_data,
                    Real a_dx,
                    Real a_time) const;

  AMRLevelMappedStabilityStrategy* clone() const;

};

#include "NamespaceFooter.H"

#endif
