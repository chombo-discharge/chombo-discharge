#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDSTABILITYSTRATEGY_H_
#define _AMRLEVELMAPPEDSTABILITYSTRATEGY_H_

#include "MOLPhysics.H"
#include "LevelGridMetrics.H"
#include "MappedLevelData.H"

#include "NamespaceHeader.H"

//! \class AMRLevelMappedStabilityStrategy
//! This class defines an interface for computing quantities needed for
//! deciding on a stable time step for high-order mapped AMR algorithms.
class AMRLevelMappedStabilityStrategy
{
  public:

  //! Base class strategy constructor. Must be called by all subclasses.
  //! \param a_stabilityFactor The factor representing the broadest height of the region
  //!                          of stability for the time integrator of interest.
  explicit AMRLevelMappedStabilityStrategy(Real a_stabilityFactor);

  /// Destructor.
  virtual ~AMRLevelMappedStabilityStrategy();

  //! This returns the maximum allowable time step (without a CFL safety factor.
  //! \param a_data The mapped grid data to be used in determining the time step.
  //! \param a_dx The grid spacing that determines the time step.
  //! \param a_time The time at which the time step is computed.
  Real maxTimeStep(const ProblemDomain& a_problem_domain,
                   const DisjointBoxLayout& a_grids,
                   const LevelGridMetrics* const a_gridMetrics,
                   const MOLPhysics& a_physics,
                   const MappedLevelData& a_data,
                   Real a_dx,
                   Real a_time) const;

  //! Override this method to determine the maximum wave speed
  //! on the computational domain.
  //! \param a_data The mapped grid data to be used in determining the wave speed.
  virtual Real maxWaveSpeed(const ProblemDomain& a_problem_domain,
                            const DisjointBoxLayout& a_grids,
                            const LevelGridMetrics* const a_gridMetrics,
                            const MOLPhysics& a_physics,
                            const MappedLevelData& a_data,
                            Real a_dx,
                            Real a_time) const = 0;

  //! Override this method to create a new AMRLevelMappedStabilityStrategy instance.
  virtual AMRLevelMappedStabilityStrategy* clone() const = 0;

  protected:

  //! The stability factor.
  Real m_stabilityFactor;

  private:

  // Default constructor, copy and assignment not allowed
  AMRLevelMappedStabilityStrategy();
  AMRLevelMappedStabilityStrategy(const AMRLevelMappedStabilityStrategy&);
  AMRLevelMappedStabilityStrategy& operator=(const AMRLevelMappedStabilityStrategy&);
};

#include "NamespaceFooter.H"

#endif
