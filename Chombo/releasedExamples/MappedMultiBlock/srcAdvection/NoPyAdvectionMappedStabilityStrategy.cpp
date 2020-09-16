#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NoPyAdvectionMappedStabilityStrategy.H"
#include "EdgeToCell.H"
// #include "FourthOrderUtil.H"
// #include "MaxScaledAcousticSpeedF_F.H"
#include "MOLUtilities.H"
// #include "SetCentersF_F.H"
// #include "BoxIterator.H"
#include "MOLAdvectionPhysics.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
NoPyAdvectionMappedStabilityStrategy::
NoPyAdvectionMappedStabilityStrategy(Real a_stabilityFactor) :
  AMRLevelMappedStabilityStrategy(a_stabilityFactor)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
NoPyAdvectionMappedStabilityStrategy::
~NoPyAdvectionMappedStabilityStrategy()
{
}

//-----------------------------------------------------------------------
Real
NoPyAdvectionMappedStabilityStrategy::
maxWaveSpeed(const ProblemDomain& a_problem_domain,
             const DisjointBoxLayout& a_grids,
             const LevelGridMetrics* const a_gridMetrics,
             const MOLPhysics& a_physics,
             const MappedLevelData& a_data,
             Real a_dx,
             Real a_time) const
{
  // Set up utilities for deconvolution.
  MOLUtilities util;
  util.define(a_problem_domain, a_dx);

  /*
  PhysAdvectMappedIBC* physIBCPtr =
    dynamic_cast<PhysAdvectMappedIBC*>(a_physics.getPhysIBC());
    
  // These lines added 29 March 2011, so that restart will work properly.
  // m_physIBCPtr->setCoordSys(m_coordSysPtr);
  physIBCPtr->setTime(a_time);

  IntVect ghostVect = m_numGhost * IntVect::Unit;
  LevelData<FArrayBox> advVel(a_grids, SpaceDim, ghostVect);
  LevelData<FluxBox> advVelFace(a_grids, SpaceDim, ghostVect);

  physIBCPtr->advVel(advVel, a_time);
  fourthOrderCellToFaceCenters(advVelFace, advVel);
  */

  Real maxVelSumLoc = 0.0;
  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  MOLAdvectionPhysics& advectionPhysics =
    dynamic_cast<MOLAdvectionPhysics&>(nonConstPhysics);
  if (a_gridMetrics != NULL)
  {
    // CH_assert(m_isAdvVelSet);
    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box box  = a_grids[dit];
      // X:  physical coordinates of faces.
      const NewFourthOrderCoordSys* coordSys = a_gridMetrics->getCoordSys(box);

      advectionPhysics.setCurrentBox(box);
      advectionPhysics.setCurrentCoordSys(coordSys);

      const FluxBox& faceAveragedVel = advectionPhysics.velFlux();

      // Box VBox = faceAveragedVel.box();

      // Velocity in computational space, given as component 'd' on face 'd'.
      // This is N^t * v.
      FluxBox faceAvgVelCSpc(box, 1);
      coordSys->computeMetricTermProductAverage(
          faceAvgVelCSpc, // Result: Velocity in CSpc
          faceAveragedVel, // Velocity in PSpc
          a_gridMetrics->m_N[dit],
          box,
          false); // 2nd order

      // Compute the maximum wave speed on cells from the face-averaged
      // velocities.
      FArrayBox maxAbsVel(box, SpaceDim);
      // Set maxWaveSpeed to sum of maxAbsVel over all dimensions.
      FArrayBox maxWaveSpeed(box, 1);
      maxWaveSpeed.setVal(0.);
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        // Set maxAbsVel(iv, dir)
        // to max{|faceAvgVelCSpc[dir](iv +/- 1/2*BASISV(dir), 0)|}.
        int faceComp = 0, cellComp = dir;
        faceAvgVelCSpc[dir].abs();
        EdgeToCellMax(faceAvgVelCSpc, faceComp,
                      maxAbsVel, cellComp, dir);
        // Divide |N^t * v| by <J> to get the true computational wave speed.
        maxAbsVel.divide(a_gridMetrics->m_J[dit], box, 0, cellComp, 1);
        maxWaveSpeed.plus(maxAbsVel, cellComp, 0);
      }

      const Real thisBoxMax = maxWaveSpeed.norm(box, 0, 0, 1);
      if (thisBoxMax > maxVelSumLoc) maxVelSumLoc = thisBoxMax;
    }
  }
  else
  {
    //**FIXME for no grid metrics (make a sum of max wave speed over Dirs.)
    const LevelData<FArrayBox>& U = a_data.getU();
    const DisjointBoxLayout& layout = U.disjointBoxLayout();
    DataIterator dit = layout.dataIterator();
    // Loop over all grids to get the maximum wave speed
    for (dit.begin(); dit.ok(); ++dit)
    {
      // Get maximum wave speed for this grid
      const Box& bx = layout[dit];
      const FArrayBox& Ufab = U[dit];
      Real speedOverBox = nonConstPhysics.getMaxWaveSpeed(Ufab, bx);
      // Compute a running maximum
      maxVelSumLoc = Max(maxVelSumLoc, speedOverBox);
    }
  }

  // Return the maximum wave speed local to this processor
  return maxVelSumLoc;
}

//-----------------------------------------------------------------------
AMRLevelMappedStabilityStrategy*
NoPyAdvectionMappedStabilityStrategy::
clone() const
{
  return new NoPyAdvectionMappedStabilityStrategy(m_stabilityFactor);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

