#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

//this one has getNearPeriodic
#include "PiecewiseLinearFillPatch.H"
#include "FourthOrderCFInterp.H"
#include "NamespaceHeader.H"

void FourthOrderCFInterp::define(/// layout at this level
                                  const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                                  /// layout at coarser level
                                  const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                                  /// number of variables
                                  const int&                a_numStates,
                                  /// problem domain on the coarser level
                                  const ProblemDomain&      a_coarseDomain,
                                  /// refinement ratio between this level and the coarser level
                                  const int&                a_refineCoarse,
                                  /// number of layers of ghost cells to fill by interpolation
                                  const int&                a_interpRadius)
{
  // Cache data
  m_numStates = a_numStates;
  m_coarseDomain = a_coarseDomain;
  m_refineCoarse = a_refineCoarse;
  m_interpRadius = a_interpRadius;

  m_layout = a_thisDisjointBoxLayout;
  m_coarseLayout = a_coarserDisjointBoxLayout;

  m_refineVect = m_refineCoarse * IntVect::Unit;

  ProblemDomain fineDomain = refine(m_coarseDomain, m_refineVect);

  // width of ghost layer of coarse cells with fine cells to fill in
  // int ghostCoarsened = ceil(m_interpRadius / m_refineCoarse);
  int ghostsCoarsened = m_interpRadius / m_refineCoarse;
  if (ghostsCoarsened * m_refineCoarse < m_interpRadius) ghostsCoarsened++;

  //no dimensions are fixed
  // FourthOrderFineInterp m_spaceInterpolator;
  m_spaceInterpolator.define(m_layout, m_numStates,
                             m_refineCoarse, fineDomain, ghostsCoarsened);

  DisjointBoxLayout m_layoutCoarsened;
  coarsen(m_layoutCoarsened, m_layout, m_refineVect);

  // LayoutData<IntVectSet> m_coarsenedGhosts;
  m_coarsenedGhosts.define(m_layoutCoarsened);

  // I copied this code segment from PiecewiseLinearFillPatch::define().
  // We find LayoutData<IntVectSet> m_coarsenedGhosts, using:
  // ProblemDomain m_coarseDomain;
  // DisjointBoxLayout m_layoutCoarsened;
  // int ghostsCoarsened;

  // a box which will determine whether a given box
  // adjoins a periodic boundary
  Box periodicTestBox(m_coarseDomain.domainBox());
  if (m_coarseDomain.isPeriodic())
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (m_coarseDomain.isPeriodic(idir))
            periodicTestBox.grow(idir,-1);
        }
    }

  bool isSorted = (a_thisDisjointBoxLayout.isSorted() && a_coarserDisjointBoxLayout.isSorted());

  DataIterator dit = m_layoutCoarsened.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& coarseBox = m_layoutCoarsened[dit];
      Box coarseBoxWithGhosts =
        grow(coarseBox, ghostsCoarsened) & m_coarseDomain;
      
      Box fine_box = m_layout[dit()];
      ProblemDomain fine_problem_domain = fineDomain;
      int interp_radius = m_interpRadius;
      bool nearPeriodic = getNearPeriodic(fine_box, fine_problem_domain,
                                          interp_radius);

      // Initially set coarsenedGhostsHere to the whole ghosted coarse box.
      IntVectSet& coarsenedGhostsHere = m_coarsenedGhosts[dit];
      coarsenedGhostsHere = IntVectSet(coarseBoxWithGhosts);

      // Iterate over boxes in coarsened fine layout, and subtract off
      // from the set of coarse cells from which the fine ghost cells
      // will be interpolated.
      LayoutIterator litOther = m_layoutCoarsened.layoutIterator();
      for (litOther.begin(); litOther.ok(); ++litOther)
        {
          const Box& other_coarsened_box
            = m_layoutCoarsened.get(litOther());
          const Box&  testBox = other_coarsened_box;

          //moving window optimizations
          if (!nearPeriodic && isSorted && (testBox.bigEnd(0) < coarseBoxWithGhosts.smallEnd(0)))
            {
              // can skip rest since  we haven't gotten
              // to something interesting
              continue;
            }

          if (!nearPeriodic && isSorted && (testBox.smallEnd(0) > coarseBoxWithGhosts.bigEnd(0)))
            {
              //can break out of loop, since we know that
              // the smallEnd of all the remaining boxes
              // are lexigraphically beyond this ghosted
              //box.
              break;
            }

          const Box& coarseOtherBox = m_layoutCoarsened[litOther];
          coarsenedGhostsHere -= coarseOtherBox;
          // also need to remove periodic images from list of cells
          // to be filled, since they will be filled through exchange
          // as well
          if (m_coarseDomain.isPeriodic()
              && !periodicTestBox.contains(coarseOtherBox)
              && !periodicTestBox.contains(coarseBox))
            {
              ShiftIterator shiftIt = m_coarseDomain.shiftIterator();
              IntVect shiftMult(m_coarseDomain.domainBox().size());
              Box shiftedBox(coarseOtherBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  coarsenedGhostsHere -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
    }

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderCFInterp::coarseFineInterp(/// interpolated solution on this level
                                           LevelData<FArrayBox>&         a_fineData,
                                           /// solution on coarser level
                                           const LevelData<FArrayBox>&   a_coarseData,
                                           /// starting coarse data component
                                           int                           a_srcComp,
                                           /// starting fine data component
                                           int                           a_dstComp,
                                           /// number of data components to interpolate
                                           int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);

  a_coarseData.copyTo(srcInterval,
                      m_spaceInterpolator.coarsenedFineData(),
                      dstInterval);

  // Interpolate to a_fineData from m_spaceInterpolator.coarsenedFineData(),
  // on the given components.
  fillInterpSpaceFromCoarsened(a_fineData, a_srcComp, a_dstComp, a_numComp);
}


//////////////////////////////////////////////////////////////////////////////
void 
FourthOrderCFInterp::
homogeneousCoarseFineInterp(/// interpolated solution on this level
                            LevelData<FArrayBox>&         a_fineData,
                            /// solution on coarser level
                            int                           a_srcComp,
                            /// starting fine data component
                            int                           a_dstComp,
                            /// number of data components to interpolate
                            int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);

  LevelData<FArrayBox>& coarseData = m_spaceInterpolator.coarsenedFineData();
  for(DataIterator dit = coarseData.dataIterator(); dit.ok(); ++dit)
    {
      coarseData[dit()].setVal(0.);
    }

  // Interpolate to a_fineData from m_spaceInterpolator.coarsenedFineData(),
  // on the given components.
  fillInterpSpaceFromCoarsened(a_fineData, a_srcComp, a_dstComp, a_numComp);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderCFInterp::fillInterpSpaceFromCoarsened(/// interpolated solution on this level
                                                        LevelData<FArrayBox>&         a_fineData,

                                                        /// starting coarse data component
                                                        int                           a_srcComp,
                                                        /// starting fine data component
                                                        int                           a_dstComp,
                                                        /// number of data components to interpolate
                                                        int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  // This should be const, but aliasLevelData doesn't let you do that.
  LevelData<FArrayBox>& coarsenedFineData =
    m_spaceInterpolator.coarsenedFineData();
  LevelData<FArrayBox> coarseCompData;
  aliasLevelData(coarseCompData, &coarsenedFineData, srcInterval);

  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);
  LevelData<FArrayBox> fineCompData;
  aliasLevelData(fineCompData, &a_fineData, dstInterval);

  DataIterator dit = m_layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Do not fill in all of fineFab:  fill ghost cells only.
      FArrayBox& fineFab = fineCompData[dit];
      const FArrayBox& coarseFab = coarseCompData[dit];
      const IntVectSet& ivs = m_coarsenedGhosts[dit];

      m_spaceInterpolator.interpOnPatch(fineFab, coarseFab, dit, ivs);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
  //  m_spaceInterpolator.interpToFine(a_fineData,
  //                                   coarsenedFineData);

  // overwrite interpolated data with valid data when present
  // a_fineData.exchange();
}



#include "NamespaceFooter.H"
