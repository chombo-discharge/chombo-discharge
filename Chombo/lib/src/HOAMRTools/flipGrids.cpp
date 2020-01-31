#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// flipGrids.cpp
//
// Dan Martin, Thurs March 5, 2009

#include "flipGrids.H"
#include "LayoutIterator.H"
#include "DataIterator.H"

#include "NamespaceHeader.H"

/// function to intersect DisjointBoxLayout with a Box, then reflect it.
/**
The getFlippedGrids function is designed to aid in the construction of
reflection boundary conditions. What it does:

For each box in a_grids,
1. intersect the grid box with a_dataBox. If a nonzero intersection, then
2. reflect that intersection Box about the reflectCoordinate in the
   reflectDirection.
3. create a DisjointBoxLayout of these boxes for which the processor
   assignments are the same as the origin boxes.
4. Create a correlation list to correlate the DataIndex's of the
   reflectedBoxes with those in the original DisjointBoxLayout.
 */
void
getFlippedGrids(DisjointBoxLayout& a_reflectedBoxes,
                Vector<Tuple<DataIndex,2> >& a_boxCorrelation,
                const DisjointBoxLayout& a_grids,
                const Box& a_dataBox,
                int a_reflectDir,
                int a_reflectCoord)
{

  Vector<Box> newBoxes;
  Vector<int> newProcAssign;

  // reset box correlation, if necessary
  a_boxCorrelation = Vector<Tuple<DataIndex, 2> >();

  LayoutIterator lit = a_grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      // intersect this grid box with the dataBox
      Box testBox = a_grids[lit];
      testBox &= a_dataBox;
      if (!testBox.isEmpty())
        {
          IntVect newLo = testBox.smallEnd();
          IntVect newHi = testBox.bigEnd();

          // flip box around reflection point
          newLo[a_reflectDir] = 2*a_reflectCoord -testBox.bigEnd(a_reflectDir);
          newHi[a_reflectDir] = 2*a_reflectCoord -testBox.smallEnd(a_reflectDir);

          Box newBox(newLo, newHi);
          newBoxes.push_back(newBox);
          newProcAssign.push_back(a_grids.procID(lit()));
        }  // end if this is a relevant box
    } // end loop over boxes in original DBL

  if (newBoxes.size() == 0)
    {
      // just set to an empty DBL and be done with it
      a_reflectedBoxes = DisjointBoxLayout();
    }
  else
    {
      // define a new DBL with the reflected boxes and corresponding
      // processor assignments
      a_reflectedBoxes.define(newBoxes, newProcAssign);
    }

  // now set up box correlations
  DataIterator dit = a_reflectedBoxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // shift it back...
      Box thisBox(a_reflectedBoxes[dit]);
      IntVect oldLo = thisBox.smallEnd();
      IntVect oldHi = thisBox.bigEnd();

      // flip box back around reflection point
      oldLo[a_reflectDir] = 2*a_reflectCoord -thisBox.bigEnd(a_reflectDir);
      oldHi[a_reflectDir] = 2*a_reflectCoord -thisBox.smallEnd(a_reflectDir);

      Box testBox(oldLo, oldHi);

      // now loop over boxes in the original DisjointBoxLayout and
      // see which one fits.
      DataIterator oldDit = a_grids.dataIterator();
      for (oldDit.begin(); oldDit.ok(); ++oldDit)
        {
          if (a_grids[oldDit].contains(testBox))
            {
              Tuple<DataIndex,2> correlationItem;
              correlationItem[0] = oldDit();
              correlationItem[1] = dit();
              a_boxCorrelation.push_back(correlationItem);

              break; // can break out of inner loop nowe
            }
        }

    } // end loop over reflected grids

  // reality check (this doesn't work in parallel, so I'm
  // commenting it out)
  //CH_assert(a_boxCorrelation.size() == a_reflectedBoxes.size());

}

#include "NamespaceFooter.H"




