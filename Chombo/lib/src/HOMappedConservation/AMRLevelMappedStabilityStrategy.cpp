#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevelMappedStabilityStrategy.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AMRLevelMappedStabilityStrategy::
AMRLevelMappedStabilityStrategy(Real a_stabilityFactor):
  m_stabilityFactor(a_stabilityFactor)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMRLevelMappedStabilityStrategy::
~AMRLevelMappedStabilityStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AMRLevelMappedStabilityStrategy::
maxTimeStep(const ProblemDomain& a_problem_domain,
            const DisjointBoxLayout& a_grids,
            const LevelGridMetrics* const a_gridMetrics,
            const MOLPhysics& a_physics,
            const MappedLevelData& a_data,
            Real a_dx,
            Real a_time) const
{
  Real maxVelSumLoc = maxWaveSpeed(a_problem_domain, a_grids, a_gridMetrics, a_physics,
                                   a_data, a_dx, a_time);
  Real maxVelSum = maxVelSumLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelSumLoc,
                &maxVelSum,
                1,
                MPI_CH_REAL,
                MPI_MAX,
                Chombo_MPI::comm);
#endif

  if (maxVelSum > 1e-15)
    return m_stabilityFactor * a_dx / maxVelSum;
  else
    return FLT_MAX;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

