#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectionGradientMappedTaggingStrategy.H"
#include "AdvectPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "LoHiCenter.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AdvectionGradientMappedTaggingStrategy::
AdvectionGradientMappedTaggingStrategy(Real a_refinementThreshold,
                                       bool a_refinementIsScaled,
                                       bool a_relativeGradient):
  AMRLevelMappedTaggingStrategy(),
  m_refineThresh(a_refinementThreshold),
  m_refinementIsScaled(a_refinementIsScaled),
  m_relativeGradient(a_relativeGradient)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AdvectionGradientMappedTaggingStrategy::
~AdvectionGradientMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
AdvectionGradientMappedTaggingStrategy::
cellsToRefine(int a_level,
              const ProblemDomain& a_problem_domain, // not used
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
  int tagIndex = 0;

  Real threshold = m_refineThresh;
  // LevelData<FArrayBox> vecMag(a_grids, 1);
  if (m_refinementIsScaled) threshold *= a_dx;

  // Compute relative gradient
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = a_grids[dit];
      const ProblemDomain& blockDomain = a_gridMetrics->blockDomain(bx);

      FArrayBox vecMagFab(bx, 1);
      // FArrayBox& vecMagFab = vecMag[dit];
      int numVecComps = SpaceDim;

      FArrayBox vecFab(bx, numVecComps);
      const FArrayBox& UFab = U[dit]; // U has at least 1 ghost layer

      Box bx1 = grow(bx, 1); // need this for gradient
      FArrayBox WFab(bx1, numW);
      nonConstPhysics.consToPrim(WFab, UFab, bx1);

      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          Box inBox = grow(bx, BASISV(idir));
          Box loBox, hiBox, centerBox, entireBox;
          int hasLo, hasHi;
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     inBox, blockDomain, idir);

          if (m_relativeGradient) // (u[i+1] - u[i-1]) / (u[i+1] + u[i-1])
            FORT_GETRELGRADF(CHF_FRA1(vecFab, idir),
                             CHF_CONST_FRA1(WFab, tagIndex),
                             CHF_CONST_INT(idir),
                             CHF_BOX(loBox),
                             CHF_CONST_INT(hasLo),
                             CHF_BOX(hiBox),
                             CHF_CONST_INT(hasHi),
                             CHF_BOX(centerBox));
          else // u[i+1] - u[i-1]
            FORT_GETGRADF(CHF_FRA1(vecFab, idir),
                          CHF_CONST_FRA1(WFab, tagIndex),
                          CHF_CONST_INT(idir),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));
        }
      FORT_MAGNITUDEF(CHF_FRA1(vecMagFab,0),
                      CHF_CONST_FRA(vecFab),
                      CHF_BOX(bx));

      // Tag where vector magnitude exceeds threshold
      BoxIterator bit(bx);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (vecMagFab(iv, 0) >= threshold)
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
AdvectionGradientMappedTaggingStrategy::
clone() const
{
  return new AdvectionGradientMappedTaggingStrategy(m_refineThresh,
                                                    m_refinementIsScaled,
                                                    m_relativeGradient);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

