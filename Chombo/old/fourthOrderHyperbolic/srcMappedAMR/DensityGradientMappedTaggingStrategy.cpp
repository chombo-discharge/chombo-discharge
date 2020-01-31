#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DensityGradientMappedTaggingStrategy.H"
#include "PolytropicPhysicsF_F.H"
#include "VorticityF_F.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
DensityGradientMappedTaggingStrategy::
DensityGradientMappedTaggingStrategy(Real a_refinementThreshold,
                                     bool a_refinementIsScaled,
                                     bool a_tagPressure,
                                     bool a_tagVorticity):
  AMRLevelMappedTaggingStrategy(),
  m_refineThresh(a_refinementThreshold),
  m_refinementIsScaled(a_refinementIsScaled),
  m_tagPressure(a_tagPressure),
  m_tagVorticity(a_tagVorticity)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
DensityGradientMappedTaggingStrategy::
~DensityGradientMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
DensityGradientMappedTaggingStrategy::
cellsToRefine(int a_level,
              const ProblemDomain& a_problem_domain,
              const DisjointBoxLayout& a_grids,
              LevelGridMetrics* a_gridMetrics,
              const MOLPhysics& a_physics,
              const MappedLevelData& a_data,
              Real a_dx,
              Real a_time) const
{
  // Create tags based on undivided gradient of density
  IntVectSet localTags;

  // Get <U>.  We need it in valid + 1 layer of (valid and invalid) ghost cells
  const LevelData<FArrayBox>& U = a_data.getU(1, 1);

  int numW = a_physics.numPrimitives();

  // Index into primitive variables.
  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  int tagIndex = m_tagPressure ? nonConstPhysics.pressureIndex()
                               : nonConstPhysics.densityIndex();

  // for vorticity, use m_refineThresh = 20.
  Real threshold = m_refineThresh;
  LevelData<FArrayBox> vecMag(a_grids, 1);
  if (m_tagVorticity)
  {
    MayDay::Error("Tagging with vorticity probably doesn't work with metrics"
        "yet");
    // also fills in vecMag
    threshold = maxVorticity(vecMag, a_gridMetrics, a_physics, a_data, a_dx);
    threshold *= m_refineThresh;
  }
  else
    {
      // petermc, 15 July 2009:  possibility for m_refinementIsScaled
      if (m_refinementIsScaled) threshold *= a_dx;
    }

  // Compute relative gradient
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = a_grids[dit];
      FArrayBox& vecMagFab = vecMag[dit];
      if (!m_tagVorticity)
        {
          // If tagging on pressure or density, vecFab is gradient.
          int numVecComps = SpaceDim;

          FArrayBox vecFab(bx, numVecComps);
          const FArrayBox& UFab = U[dit];

          Box bx1 = grow(bx, 1); // need this for gradient
          FArrayBox WFab(bx1, numW);
          nonConstPhysics.consToPrim(WFab, UFab, bx1);

          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              const Box bxCenter = bx & grow(a_problem_domain,-BASISV(idir));

              const Box bxLo     = bx & adjCellLo(bxCenter,idir);
              const int hasLo = ! bxLo.isEmpty();

              const Box bxHi     = bx & adjCellHi(bxCenter,idir);
              const int hasHi = ! bxHi.isEmpty();

              FORT_GETRELGRADF(CHF_FRA1(vecFab, idir),
                               CHF_CONST_FRA1(WFab, tagIndex),
                               CHF_CONST_INT(idir),
                               CHF_BOX(bxLo),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(bxHi),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(bxCenter));
            }
          FORT_MAGNITUDEF(CHF_FRA1(vecMagFab,0),
                          CHF_CONST_FRA(vecFab),
                          CHF_BOX(bx));
        }

      // Tag where vector magnitude exceeds threshold
      BoxIterator bit(bx);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          if (vecMagFab(iv) >= threshold)
            {
              localTags |= iv;
            }
        }
    }

  // End application-dependent code - PC.
  return localTags;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
DensityGradientMappedTaggingStrategy::
maxVorticity(LevelData<FArrayBox>& a_vecMag,
             LevelGridMetrics* a_gridMetrics,
             const MOLPhysics& a_physics,
             const MappedLevelData& a_data,
             Real a_dx) const
{
  // Get <U>.  We need it in valid + 1 layer of (valid and invalid) ghost cells
  const LevelData<FArrayBox>& U = a_data.getU(1, 1);

  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  int numW = nonConstPhysics.numPrimitives();
  const Interval& velInt = nonConstPhysics.velocityInterval();
  int numVecComps = SpaceDim;
  if (SpaceDim == 2) numVecComps = 1;
  Real vorticityMaxLocal = 0.;
  const DisjointBoxLayout& levelDomain = a_gridMetrics->getBoxes();
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = levelDomain[dit];
      FArrayBox vecFab(bx, numVecComps);
      const FArrayBox& UFab = U[dit];

      Box bx1 = grow(bx, 1); // need this for gradient
      FArrayBox WFab(bx1, numW);
      FArrayBox velFab(velInt, WFab);
      nonConstPhysics.consToPrim(WFab, UFab, bx1);
      // For tagging purposes, we don't need highly accurate
      // estimate of vorticity.
      if (SpaceDim == 3)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              FORT_COMPUTEVORT(CHF_FRA1(vecFab, idir),
                               CHF_CONST_FRA(velFab),
                               CHF_BOX(bx),
                               CHF_CONST_REAL(a_dx),
                               CHF_CONST_INT(idir));
            }
        }
      else if (SpaceDim == 2)
        {
          int idir = 0; // ignored by FORT_COMPUTEVORT
          FORT_COMPUTEVORT(CHF_FRA1(vecFab, idir),
                           CHF_CONST_FRA(velFab),
                           CHF_BOX(bx),
                           CHF_CONST_REAL(a_dx),
                           CHF_CONST_INT(idir));
        }
      else
        {
          MayDay::Error("Vorticity requires 2D or 3D");
        }

      FArrayBox& vecMagFab = a_vecMag[dit];
      FORT_MAGNITUDEF(CHF_FRA1(vecMagFab, 0),
                      CHF_CONST_FRA(vecFab),
                      CHF_BOX(bx));
      Real fabMax = vecMagFab.max();
      if (fabMax > vorticityMaxLocal) vorticityMaxLocal = fabMax;
    }
  Real vorticityMax = vorticityMaxLocal;
#ifdef CH_MPI
  MPI_Allreduce(&vorticityMaxLocal, &vorticityMax, 1,
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return vorticityMax;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMRLevelMappedTaggingStrategy*
DensityGradientMappedTaggingStrategy::
clone() const
{
  return new DensityGradientMappedTaggingStrategy(m_refineThresh,
               m_refinementIsScaled, m_tagPressure, m_tagVorticity);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

