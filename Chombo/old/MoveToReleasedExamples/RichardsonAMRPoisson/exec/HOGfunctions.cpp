#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"

#include "localfunctions.H"


void createHOG(Vector<LDFB*>& hog,
               const int numLevels,
               const Vector<DBL>& DBLVector)
{
 CH_assert(hog.size() >= numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    hog[ilev] = new LDFB(DBLVector[ilev], 1, unitVector);
  }
}

void multiplyHOGbyConstant(Vector<LDFB*>& hog,
                           const int numLevels,
                           const Real k1)
{
 CH_assert(hog.size() >= numLevels);
  for (int ilev=0; ilev<numLevels; ++ilev)
  {
    LDFB& levhog = *hog[ilev];
    for (DataIterator dit = levhog.dataIterator(); dit.ok(); ++dit)
    {
      levhog[dit()] *= k1;
    }
  }
}

void plusEqualHOG(Vector<LDFB*>& hog1,
                  const int numLevels,
                  const Vector<LDFB*>& hog2)
{
 CH_assert(hog1.size() >= numLevels);
 CH_assert(hog2.size() >= numLevels);
  for (int ilev=0; ilev<numLevels; ++ilev)
  {
    LDFB& levhog1 = *hog1[ilev];
    const LDFB& levhog2 = *hog2[ilev];
    for (DataIterator dit = levhog1.dataIterator(); dit.ok(); ++dit)
    {
      levhog1[dit()] += levhog2[dit()];
    }
  }
}



void copyFirstHOGintoSecond(const Vector<LDFB*>& first,
                            Vector<LDFB*>& second,
                            const int numLevels)
{
  for (int ilev=0; ilev<numLevels; ++ilev)
  {
    LDFB& lsecond = *second[ilev];
    const LDFB& lfirst = *first[ilev];
    lfirst.copyTo(lfirst.interval(), lsecond, lsecond.interval());
  }
}

void subtractTwoHOGs(const Vector<LDFB*>& grid1,
                     const Vector<LDFB*>& grid2,
                     const int numLevels,
                     Vector<LDFB*>& difference)
{

 CH_assert(grid1.size() >= numLevels);
 CH_assert(grid1.size() >= numLevels);
 CH_assert(difference.size() >= numLevels);

  copyFirstHOGintoSecond(grid1, difference, numLevels);
  //for (int ilev=0; ilev<numLevels; ++ilev)
  //{
  //LDFB& levdiff = *difference[ilev];
  //const LDFB& levgrid1 = *grid1[ilev];
  ////difference[ilev] = grid1[ilev];
  //levgrid1.copyTo(levgrid1.interval(), levdiff, levdiff.interval());
  //}

  for (int ilev=0; ilev<numLevels; ++ilev)
  {
    LDFB& levdifference = *difference[ilev];
    const LDFB& lev_grid2 = *grid2[ilev];
    for (DataIterator dit = levdifference.dataIterator(); dit.ok(); ++dit)
    {
      levdifference[dit()] -= lev_grid2[dit()];
    }
  }
}


void averageDownHOG(Vector<LDFB*>& grid,
                    const Vector<int>& refinementRatioVector,
                    const int numLevels)
{

  if (numLevels == 0) return;
 CH_assert(grid.size() >= numLevels);
 CH_assert(refinementRatioVector.size() >= numLevels);

  // Start at the second to the finest level:  numLevels-2
  for (int i=numLevels-2; i>=0; i--)
  {

    const DBL currentDBL = (*grid[i]).getBoxes();
    const DBL finerDBL = (*grid[i+1]).getBoxes();

    LDFB& currentLDFB     = *grid[i];
    const LDFB& finerLDFB = *grid[i+1];

    LDFB modifiedFinerLDFB(finerDBL, 1, unitVector);

    // overwrite the values of the modifiedFinerLDFB with those
    // from higher levels contained in finerLDFB.
    //modifiedFinerLDFB = finerLDFB;
    finerLDFB.copyTo(finerLDFB.interval(), modifiedFinerLDFB,
                    modifiedFinerLDFB.interval());

    // multiply modifiedFinerLDFB by Nref^2
    for (DataIterator dit = modifiedFinerLDFB.dataIterator(); dit.ok(); ++dit)
    {
      modifiedFinerLDFB[dit()] *= refinementRatioVector[i]*refinementRatioVector[i];
    }

    // Average down the now modified Finer LDFB onto the current LDFB

    CoarseAverage avgDown(finerDBL, 1, refinementRatioVector[i]);

    // take level in 2nd argument and average down into 1st argument
    //  so 1st argument should be coarser.
    //  1st argument is changed.
    avgDown.averageToCoarse(currentLDFB, modifiedFinerLDFB);

  }

}


void compareTwoHOGs(const Vector<LDFB*>& grid1, const Vector<LDFB*>& grid2,
                    const int numLevels, const bool returnRawDiff,
                    Vector<LDFB*>& absDifference)
{

  //printf(" grid1 size = %d   grid2 size = %d\n", grid1.size(), grid2.size());
 CH_assert(grid1.size() >= numLevels);
 CH_assert(grid1.size() >= numLevels);
 CH_assert(absDifference.size() >= numLevels);

  // Subtract grid2 from grid1 and store into diff
  Vector<LDFB*> diff(numLevels, NULL);
  subtractTwoHOGs(grid1, grid2, numLevels, diff);
  //prettyPrintLDFB(diff[0], " zero?");

  // now find the largest absolute value in matrix
  copyFirstHOGintoSecond(absDifference, diff, numLevels);

  if (returnRawDiff) return;
  // no need to do the rest

  Vector<Real> maxError(numLevels);
  Vector<IntVect> indexOfMaxError(numLevels);

  for (int i=0; i<numLevels; ++i)
  {
    LDFB& levDiff = *(absDifference[i]);
    const DBL dbl = levDiff.disjointBoxLayout();
    maxError[i] = 0;

    for (DataIterator dit = levDiff.dataIterator(); dit.ok(); ++dit)
    {
      levDiff[dit()].abs();

      //printf(" max=%16.7e at (%4d,%4d)  min=%16.7e at (%4d,%4d)  norm=%16.7e n=%ld\n",
      //   levDiff[dit()].max(dbl[dit()]),
      //   (levDiff[dit()].maxIndex(dbl[dit()]))[0],
      //   (levDiff[dit()].maxIndex(dbl[dit()]))[1],
      //   levDiff[dit()].min(dbl[dit()]),
      //   (levDiff[dit()].minIndex(dbl[dit()]))[0],
      //   (levDiff[dit()].minIndex(dbl[dit()]))[1],
      //   levDiff[dit()].norm(dbl[dit()]),
      //   dbl[dit()].numPts());

      if (levDiff[dit()].max(dbl[dit()]) > maxError[i])
      {
        maxError[i] = levDiff[dit()].max(dbl[dit()]);
        indexOfMaxError[i] =  levDiff[dit()].maxIndex(dbl[dit()]);
      }

    }
  }

  //printf(" size=%10.0f  h=%16.7e  max_error=%16.7e  at (%d, %d)\n",
  // 1.0/DxVector[0], DxVector[0], maxError,
  // indexOfMaxError[0], indexOfMaxError[1]);
  for (int i=0; i<numLevels; i++)
  {
    printf("  level %d max_error=%16.7e  at (%d, %d)\n",
           i, maxError[i],
           indexOfMaxError[i][0], indexOfMaxError[i][1]);
  }

}

