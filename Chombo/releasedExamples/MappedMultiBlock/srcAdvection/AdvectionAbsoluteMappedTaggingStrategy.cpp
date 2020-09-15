#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectionAbsoluteMappedTaggingStrategy.H"
#include "AdvectPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "LoHiCenter.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AdvectionAbsoluteMappedTaggingStrategy::
AdvectionAbsoluteMappedTaggingStrategy(Real a_refinementThreshold) :
  AMRLevelMappedTaggingStrategy(),
  m_refineThresh(a_refinementThreshold)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AdvectionAbsoluteMappedTaggingStrategy::
~AdvectionAbsoluteMappedTaggingStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
IntVectSet
AdvectionAbsoluteMappedTaggingStrategy::
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

  // Get <U>.  We do not need it on any ghost cells.
  const LevelData<FArrayBox>& U = a_data.getU();

  int numW = a_physics.numPrimitives();

  // Index into primitive variables.
  MOLPhysics& nonConstPhysics = const_cast<MOLPhysics&>(a_physics); // FIXME: Stupid
  int tagIndex = 0;

  // Compute relative gradient
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bx = a_grids[dit];

      const FArrayBox& UFab = U[dit];
      FArrayBox WFab(bx, numW);
      nonConstPhysics.consToPrim(WFab, UFab, bx);

      // Tag where W[tagIndex] exceeds threshold
      BoxIterator bit(bx);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (abs(WFab(iv, tagIndex)) >= m_refineThresh)
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
AdvectionAbsoluteMappedTaggingStrategy::
clone() const
{
  return new AdvectionAbsoluteMappedTaggingStrategy(m_refineThresh);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

