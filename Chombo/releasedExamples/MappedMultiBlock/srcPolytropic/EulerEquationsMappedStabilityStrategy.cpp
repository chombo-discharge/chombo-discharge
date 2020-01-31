#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EulerEquationsMappedStabilityStrategy.H"
#include "CellToEdge.H"
#include "PolytropicPhysicsF_F.H"
#include "EdgeToCellF_F.H"
#include "FourthOrderUtil.H"
#include "MaxScaledAcousticSpeedF_F.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EulerEquationsMappedStabilityStrategy::
EulerEquationsMappedStabilityStrategy(Real a_stabilityFactor):
  AMRLevelMappedStabilityStrategy(a_stabilityFactor)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EulerEquationsMappedStabilityStrategy::
~EulerEquationsMappedStabilityStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
EulerEquationsMappedStabilityStrategy::
maxWaveSpeed(const ProblemDomain& a_problem_domain,
             const DisjointBoxLayout& a_grids,
             const LevelGridMetrics* const a_gridMetrics,
             const MOLPhysics& a_physics,
             const MappedLevelData& a_data,
             Real a_dx,
             Real a_time) const
{
  Real maxVelSumLoc = 0.0;
  int numPrim = a_physics.numPrimitives();
  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  if (a_gridMetrics != NULL)
  {
    // Need 1 layer of ghost cells
    const LevelData<FArrayBox>& cellAvgULev = a_data.getU(1, 1);

    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      //**FIXME Better to use a 1-sided op at domain boundaries for
      //**cellToEdge instead of extrapolating
      const Box box  = a_grids[dit];
      const Box box1 = grow(box, 1);
      const Box box1Dom = box1 & a_problem_domain;
      const FArrayBox& cellAvgU = cellAvgULev[dit];
      // Convert to primitive state
      FArrayBox cellAvgW(box1, numPrim);
      nonConstPhysics.consToPrim(cellAvgW, cellAvgU, box1Dom);
      // Second-order extrapolation at domain boundaries
      secondOrderCellExtrapAtDomainBdry(cellAvgW,
          box,
          a_problem_domain);
      // Cell-to-edge
      FluxBox faceAvgW(box, numPrim);
      CellToEdge(cellAvgW, faceAvgW);

      // Velocity in computational space, given as component 'd' on face 'd'
      FluxBox faceAvgVelCSpc(box, 1);
      Interval velocityInt = nonConstPhysics.velocityInterval();
      CH_assert(velocityInt.begin() <= velocityInt.end());
      // Interval velIntv(nonConstPhysics.velocityInterval().begin(),
      //                  nonConstPhysics.velocityInterval().begin());
      // a_gridMetrics->getCoordSys(box)
      //   ->computeMetricTermProductAverage(
      //       faceAvgVelCSpc,                 // Result: Velocity in CSpc
      //       faceAvgW,                       // Velocity in PSpc
      //       a_gridMetrics->m_N[dit],
      //       SpaceDim,
      //       faceAvgW,                       // Unused for 2nd order
      //       box,
      //       false,                          // 2nd order
      //       Interval(0, 0),                 // Interval for faceAvgVelCSpc
      //       velIntv);                       // Interval of Vel. in PSpc

      // petermc, 5 Sep 2012:
      // Ugh, there is no consistent FSpaceStride we can use,
      // so we must have varFStart == 0.
      // This means a_varIntervalF == (0, 0).
      // We cannot alias a FluxBox, so we must copy!
      FluxBox velAvgW(box, SpaceDim);
      velAvgW.copy(faceAvgW, velocityInt.begin(), 0, SpaceDim);
      a_gridMetrics->getCoordSys(box)
        ->computeMetricTermProductAverage(
            faceAvgVelCSpc,                 // Result: Velocity in CSpc
            velAvgW,                        // Velocity in PSpc
            a_gridMetrics->m_N[dit],
            SpaceDim,
            velAvgW,                        // Unused for 2nd order
            box,
            false,                          // 2nd order
            Interval(0, 0),                 // Interval for faceAvgVelCSpc
            Interval(0, 0));                // Interval of Vel. in PSpc

      // Face ratios
      FluxBox faceNMag(box, 1);
      a_gridMetrics->getCoordSys(box)
        ->magnitudeN(faceNMag, a_gridMetrics->m_N[dit], box);

      // Reusing cellAvgW as temp space for sum of directional wavespeeds
      cellAvgW.setVal(0., 0);
      for (int dir = 0; dir != SpaceDim; ++dir)
      {
        Box faceBox(box);
        faceBox.surroundingNodes(dir);
        FArrayBox faceSound(faceBox, 1);
        // from PolytropicPhysicsF.ChF, the speed of sound
        FORT_SOUNDSPEEDPRIMF(CHF_FRA1(faceSound, 0),
            CHF_CONST_FRA(faceAvgW[dir]),
            CHF_BOX(faceBox));

        // from MaxScaledAcousticSpeedF.ChF, scale the speed of sound
        // by the face area ratios and add it to |fluid speed| to get the
        // maximum wavespeed.
        FORT_MAXSCALEDACOUSTICSPEED(CHF_FRA1(faceAvgVelCSpc[dir], 0),
            CHF_CONST_FRA1(faceSound, 0),
            CHF_CONST_FRA1(faceNMag[dir], 0),
            CHF_BOX(faceBox));

        // Add an average of the max wave speed from the faces to the
        // cell.
        FORT_EDGETOINCREMENTCELL(CHF_CONST_FRA1(faceAvgVelCSpc[dir], 0),
            CHF_FRA1(cellAvgW, 0),
            CHF_BOX(box),
            CHF_CONST_INT(dir));
      }

      // Divide by <J> to get the true computational wave speed
      cellAvgW.divide(a_gridMetrics->m_J[dit], box, 0, 0, 1);

      const Real thisBoxMax = cellAvgW.norm(box, 0, 0, 1);
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

//-----------------------------------------------------------------------
AMRLevelMappedStabilityStrategy*
EulerEquationsMappedStabilityStrategy::
clone() const
{
  return new EulerEquationsMappedStabilityStrategy(m_stabilityFactor);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

