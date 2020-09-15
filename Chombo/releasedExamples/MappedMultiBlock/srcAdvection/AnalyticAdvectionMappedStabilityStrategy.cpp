#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AnalyticAdvectionMappedStabilityStrategy.H"
#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
#include "MOLUtilities.H"
#include "SetCentersF_F.H"
#include "BoxIterator.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AnalyticAdvectionMappedStabilityStrategy::
AnalyticAdvectionMappedStabilityStrategy(Real a_stabilityFactor,
                                         RefCountedPtr<VectorFunction> a_velocity):
  AMRLevelMappedStabilityStrategy(a_stabilityFactor),
  m_velocity(a_velocity),
  m_streamFunc()
{
  CH_assert(!a_velocity.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticAdvectionMappedStabilityStrategy::
AnalyticAdvectionMappedStabilityStrategy(Real a_stabilityFactor,
                                         RefCountedPtr<ScalarFunction> a_streamFunc):
  AMRLevelMappedStabilityStrategy(a_stabilityFactor),
  m_velocity(),
  m_streamFunc(a_streamFunc)
{
  CH_assert(!a_streamFunc.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticAdvectionMappedStabilityStrategy::
~AnalyticAdvectionMappedStabilityStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AnalyticAdvectionMappedStabilityStrategy::
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

  Real maxVelSumLoc = 0.0;
  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  if (a_gridMetrics != NULL)
  {
    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box box  = a_grids[dit];

      // Compute the face-averaged normal velocity from our velocity function.
      // We do this by computing face-centered values of V and then
      // deconvolving to obtain face-averaged values.
      FluxBox faceCenteredVel(box, 1),
              faceAveragedVel(box, 1);
      faceCenteredVel.setVal(0.0);
      faceAveragedVel.setVal(0.0);
      Box VBox = faceCenteredVel.box();

      // Xi:  cartesian coordinates
      FArrayBox XiF(VBox, SpaceDim);
      for (int d = 0; d < SpaceDim; ++d)
      {
        FORT_SETFACECENTERS(CHF_FRA(XiF),
            CHF_CONST_REAL(a_dx),
            CHF_CONST_INT(d),
            CHF_BOX(box));
      }

      // X:  physical coordinates of faces.
      const NewFourthOrderCoordSys* coordSys = a_gridMetrics->getCoordSys(box);
      int dimension = coordSys->dimension();
      FArrayBox XF(VBox, dimension);
      coordSys->realCoord(XF, XiF, VBox);

      // Compute face-centered velocities.
      if (!m_velocity.isNull())
      {
        for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          IntVect i = bit();
          RealVect X, Xii;
          for (int d = 0; d < SpaceDim; ++d)
          {
            X[d] = XF(i, d);
            Xii[d] = XiF(i, d);
          }
          RealVect U = (*m_velocity)(X, a_time);
          // Map from real to computational space.
          RealVect V = coordSys->vectorTransformationRealToMapped(Xii).transform(U); // Mapped velocity at Xii.
          for (int d = 0; d < SpaceDim; ++d)
            faceCenteredVel[d](i, 0) = V[d];
        }
      }
      else
      {
        // FIXME:
        MayDay::Error("Implement stream function velocity calculation!");
      }

      // Deconvolve to get face-averaged velocities.
      faceAveragedVel.copy(faceCenteredVel);
      util.deconvolveFace(faceAveragedVel, faceCenteredVel);

      // Velocity in computational space, given as component 'd' on face 'd'
      FluxBox faceAvgVelCSpc(VBox, 1);
      coordSys->computeMetricTermProductAverage(
          faceAvgVelCSpc, // Result: Velocity in CSpc
          faceAveragedVel, // Velocity in PSpc
          a_gridMetrics->m_N[dit],
          VBox,
          false); // 2nd order

      // Compute the maximum wave speed on cells from the face-averaged
      // velocities.
      FArrayBox maxWaveSpeed(VBox, 1);
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        int faceComp = 0, cellComp = 0;
        faceAvgVelCSpc[dir].abs();
        EdgeToCellMax(faceAvgVelCSpc, faceComp, maxWaveSpeed, cellComp, dir);

        // Divide by <J> to get the true computational wave speed.
        maxWaveSpeed.divide(a_gridMetrics->m_J[dit], VBox, 0, 0, 1);
      }

      const Real thisBoxMax = maxWaveSpeed.norm(VBox, 0, 0, 1);
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
AnalyticAdvectionMappedStabilityStrategy::
clone() const
{
  if (!m_velocity.isNull())
    return new AnalyticAdvectionMappedStabilityStrategy(m_stabilityFactor, m_velocity);
  else
    return new AnalyticAdvectionMappedStabilityStrategy(m_stabilityFactor, m_streamFunc);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

