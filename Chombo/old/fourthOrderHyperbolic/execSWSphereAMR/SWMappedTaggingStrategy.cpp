#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SWMappedTaggingStrategy.H"
#include "AdvectPhysicsF_F.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
SWMappedTaggingStrategy::
SWMappedTaggingStrategy(Real a_refinementThreshold,
                        bool a_refinementIsScaled):
  AMRLevelMappedTaggingStrategy(),
  m_refineThresh(a_refinementThreshold),
  m_refinementIsScaled(a_refinementIsScaled)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SWMappedTaggingStrategy::
~SWMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
SWMappedTaggingStrategy::
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
  int tagIndex = 0;

  Real threshold = m_refineThresh;
  LevelData<FArrayBox> vecMag(a_grids, 1);
  if (m_refinementIsScaled) threshold *= a_dx;

  // Compute relative gradient
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& bx = a_grids[dit];
    FArrayBox& vecMagFab = vecMag[dit];
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
AMRLevelMappedTaggingStrategy*
SWMappedTaggingStrategy::
clone() const
{
  return new SWMappedTaggingStrategy(m_refineThresh,
                                     m_refinementIsScaled);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

