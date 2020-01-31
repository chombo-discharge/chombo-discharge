#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelSWSourceTerm.H"
#include "ShallowWaterPhysicsF_F.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

LevelSWSourceTerm::LevelSWSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

LevelSWSourceTerm::~LevelSWSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void LevelSWSourceTerm::define(MultiBlockCoordSys* const a_coordSysPtr,
                               const MOLPhysics* const a_molPhysics,
                               const DisjointBoxLayout& a_grids)
{
  LevelSourceTerm::define(a_coordSysPtr, a_molPhysics, a_grids);
}

//////////////////////////////////////////////////////////////////////////////

LevelSourceTerm* LevelSWSourceTerm::new_sourceTerm() const
{
  LevelSourceTerm* retval = static_cast<LevelSourceTerm*>
    (new LevelSWSourceTerm());
  return retval;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelSWSourceTerm::addSourceTerm(LevelData<FArrayBox>&   a_rhs,
                                 LevelData<FArrayBox>&   a_U)
{
  // a_rhs += <J*S(a_U)>
  // a_U contains <U> and should not be modified.
  int numPrim = m_molPhysics->numPrimitives();
  int numFluxes = m_molPhysics->numFluxes();
  int numConserved = m_molPhysics->numConserved();
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rhsFab = a_rhs[dit];
      const FArrayBox& UavgFab = a_U[dit];

      const Box& baseBox = m_grids[dit];
      Box bx1 = grow(baseBox, 1);
      Box bx2 = grow(baseBox, 2);

      // Convert cell-averaged UavgFab to cell-centered UcenFab.
      FArrayBox UcenFab(bx2, numConserved);
      UcenFab.copy(UavgFab);
      fourthOrderAverageCell(UcenFab, -1); // fills in UcenFab on bx1

      int blockNumber = m_coordSysPtr->whichBlock(baseBox);
      const NewCoordSys* thisCoordSysPtr =
        m_coordSysPtr->getCoordSys(blockNumber);

      FArrayBox WcenFab(bx1, numPrim);
      ((MOLPhysics*) m_molPhysics)->consToPrim(WcenFab, UcenFab, bx1);

      FArrayBox xiFab(bx1, SpaceDim);
      thisCoordSysPtr->getCenterMappedCoordinates(xiFab, bx1);

      // Need sourceFab on bx1 only.
      // SWADDMETRICSOURCE and SWADDCORIOLISSOURCE do not need neighbors.
      FArrayBox sourceFab(bx1, numFluxes);
      sourceFab.setVal(0.);

      FORT_SWADDMETRICSOURCE(CHF_FRA(sourceFab),
                             CHF_CONST_FRA(xiFab),
                             CHF_CONST_FRA(WcenFab),
                             CHF_BOX(bx1));
      FORT_SWADDCORIOLISSOURCE(CHF_FRA(sourceFab),
                               CHF_CONST_FRA(xiFab),
                               CHF_CONST_FRA(UcenFab),
                               CHF_CONST_INT(blockNumber),
                               CHF_BOX(bx1));

      // sourceFab *= JFab
      FArrayBox JFab(bx1, 1);
      thisCoordSysPtr->pointwiseJ(JFab, xiFab, bx1);
      for (int comp = 0; comp < numFluxes; comp++)
        {
          // sourceFab[comp] *= JFab[0]
          sourceFab.mult(JFab, 0, comp);
        }

      // Convert sourceFab from cell-centered to cell-averaged.
      fourthOrderAverageCell(sourceFab);

      // rhsFab += sourceFab on baseBox, on numFluxes comps starting with 0
      rhsFab.plus(sourceFab, baseBox, 0, 0, numFluxes);
      // dummy statement in order to get around gdb bug
      int dummy_unused; dummy_unused = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
