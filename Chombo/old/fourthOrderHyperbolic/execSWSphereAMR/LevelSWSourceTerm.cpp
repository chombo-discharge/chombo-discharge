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
  // FIXME:  FOR NOW, DO NO AVERAGING.
  int numPrim = m_molPhysics->numPrimitives();
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rhsFab = a_rhs[dit];
      const FArrayBox& UFab = a_U[dit];
      const Box& baseBox = m_grids[dit];

      int blockNumber = m_coordSysPtr->whichBlock(baseBox);
      const NewCoordSys* thisCoordSysPtr =
        m_coordSysPtr->getCoordSys(blockNumber);

      FArrayBox WFab(baseBox, numPrim);
      ((MOLPhysics*) m_molPhysics)->consToPrim(WFab, UFab, baseBox);

      FArrayBox xiFab(baseBox, SpaceDim);
      thisCoordSysPtr->getCenterMappedCoordinates(xiFab, baseBox);

      FORT_SWADDMETRICSOURCE(CHF_FRA(rhsFab),
                             CHF_CONST_FRA(xiFab),
                             CHF_CONST_FRA(WFab),
                             CHF_BOX(baseBox));
      FORT_SWADDCORIOLISSOURCE(CHF_FRA(rhsFab),
                               CHF_CONST_FRA(xiFab),
                               CHF_CONST_FRA(UFab),
                               CHF_CONST_INT(blockNumber),
                               CHF_BOX(baseBox));

      // rhsFab *= JavgFab
      // Do I have only cell-averaged J and no cell-centered J?
      FArrayBox JavgFab(baseBox, 1);
      thisCoordSysPtr->getAvgJ(JavgFab, baseBox);
      for (int comp = 0; comp < numPrim; comp++)
        {
          // rhsFab[comp] *= JavgFab[0]
          rhsFab.mult(JavgFab, 0, comp);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceFooter.H"
