#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectionCosineBellMappedTaggingStrategy.H"
#include "AdvectPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "LoHiCenter.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AdvectionCosineBellMappedTaggingStrategy::
AdvectionCosineBellMappedTaggingStrategy(
                                         Real a_ambientDensity,
                                         Real a_deltaDensity,
                                         Real a_size,
                                         Real a_advectVelocity,
                                         Real a_advectAngle,
                                         RealVect a_center,
                                         Real a_refinementThreshold) :
  AMRLevelMappedTaggingStrategy(),
  m_ambientDensity(a_ambientDensity),
  m_deltaDensity(a_deltaDensity),
  m_size(a_size),
  m_advectVelocity(a_advectVelocity),
  m_advectAngle(a_advectAngle),
  m_center(a_center),
  m_refineThresh(a_refinementThreshold)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AdvectionCosineBellMappedTaggingStrategy::
~AdvectionCosineBellMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
AdvectionCosineBellMappedTaggingStrategy::
cellsToRefine(int a_level,
              const ProblemDomain& a_problem_domain, // not used
              const DisjointBoxLayout& a_grids,
              LevelGridMetrics* a_gridMetrics,
              const MOLPhysics& a_physics,
              const MappedLevelData& a_data,
              Real a_dx,
              Real a_time) const
{
  // This is copied from
  // Chombo/example/fourthOrderMappedGrids/execAdvection/AdvectCosineBell.cpp

  // Create tags based on undivided gradient of density
  IntVectSet localTags;

  // Get <U>.  We do not need it on any ghost cells.
  // const LevelData<FArrayBox>& U = a_data.getU();

  // int numW = a_physics.numPrimitives();

  // Index into primitive variables.
  // MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  // int tagIndex = 0;

  const MultiBlockCoordSys& coordSys = a_gridMetrics->getCoordSys();

  Real InitialTheta = a_time * m_advectVelocity;

  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box baseBox = a_grids[dit];

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(
          coordSys.getCoordSys(baseBox));

      CH_assert(coordSysBlockPtr);

      // Need one layer of ghost cells in order to do averaging.
      Box bx1 = grow(baseBox, 1);

      // Xi: mapped space coordinates
      FArrayBox XiFab(bx1, SpaceDim);
      coordSysBlockPtr->getCellMappedCoordinates(XiFab, bx1);

      // lonlat: longitude and latitude
      FArrayBox lonlatFab(bx1, SpaceDim);

      // FArrayBox WFab(baseBox, 1);

      // Tag where W[tagIndex] exceeds threshold
      BoxIterator bit(bx1);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          // Get equiangular coordinates
          RealVect xi;
          xi[0] = XiFab(iv,0);
          xi[1] = XiFab(iv,1);

          Real xyz[3];
          coordSysBlockPtr->pointTransformEquiangularToCartesian(xi, xyz);

          // Rotated sample points
          Real xyzR[3];

          xyzR[0] =
              xyz[0] * (sin(m_advectAngle) * sin(m_advectAngle)
                + cos(m_advectAngle) * cos(m_advectAngle) * cos(InitialTheta))
            + xyz[1] * cos(m_advectAngle) * sin(InitialTheta)
            - xyz[2] * sin(m_advectAngle) * cos(m_advectAngle)
                * (1.0 - cos(InitialTheta));

          xyzR[1] =
            - xyz[0] * cos(m_advectAngle) * sin(InitialTheta)
            + xyz[1] * cos(InitialTheta)
            - xyz[2] * sin(m_advectAngle) * sin(InitialTheta);

          xyzR[2] =
            - xyz[0] * sin(m_advectAngle) * cos(m_advectAngle)
                * (1.0 - cos(InitialTheta))
            + xyz[1] * sin(m_advectAngle) * sin(InitialTheta)
            + xyz[2] * (cos(m_advectAngle) * cos(m_advectAngle)
                + sin(m_advectAngle) * sin(m_advectAngle) * cos(InitialTheta));

          // Sample points in RLL coordinates:  longitude, then latitude
          lonlatFab(iv, 0) = atan2(xyzR[1], xyzR[0]); // - M_PI / 2.0;
          lonlatFab(iv, 1) = asin(xyzR[2]);
        }
      FArrayBox WFab(bx1, 1);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          Real dLon = lonlatFab(iv, 0);
          Real dLat = lonlatFab(iv, 1);

          // Great circle distance
          Real dR1 = cos(dLat) * sin(dLon - m_center[1]);
          Real dR2 = cos(m_center[0]) * sin(dLat)
                   - sin(m_center[0]) * cos(dLat) * cos(dLon - m_center[1]);
          Real dR3 = sin(m_center[0]) * sin(dLat)
                   + cos(m_center[0]) * cos(dLat) * cos(dLon - m_center[1]);

          Real dR = atan2(sqrt(dR1 * dR1 + dR2 * dR2), dR3);

          WFab(iv, 0) = m_ambientDensity;
          if (dR <= m_size)
            {
              // fac(dR) is between 0 and 1; fac(0) = 0; fac(m_size) = 1.
              Real fac = (1.0 + cos(M_PI * dR / m_size)) / 2.0;
              WFab(iv, 0) += m_deltaDensity * fac*fac*fac;
              // UFab(iv,0) += m_deltaDensity * fac*fac;
            }
        }

      // Change WFab from cell-centered to cell-averaged.
      ProblemDomain domain = a_gridMetrics->blockDomain(baseBox, 1);
      fourthOrderAverageCell(WFab, domain, baseBox);

      BoxIterator bitBase(baseBox);
      for (bitBase.begin(); bitBase.ok(); ++bitBase)
        {
          IntVect iv = bitBase();
          if (abs(WFab(iv, 0)) >= m_refineThresh)
            {
              localTags |= iv;
            }
        }
      // pout() << "localTags size " << localTags.numPts() << endl;
    }

  // End application-dependent code - PC.
  // pout() << "localTags size " << localTags.numPts() << endl;
  return localTags;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMRLevelMappedTaggingStrategy*
AdvectionCosineBellMappedTaggingStrategy::
clone() const
{
  return new AdvectionCosineBellMappedTaggingStrategy(
                                                      m_ambientDensity,
                                                      m_deltaDensity,
                                                      m_size,
                                                      m_advectVelocity,
                                                      m_advectAngle,
                                                      m_center,
                                                      m_refineThresh);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

