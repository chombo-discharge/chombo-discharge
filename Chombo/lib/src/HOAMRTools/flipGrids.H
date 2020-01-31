#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// flipGrids.H
//
// Dan Martin, Thurs March 5, 2009

#ifndef _FLIPGRIDS_H_
#define _FLIPGRIDS_H_

#include "DisjointBoxLayout.H"
#include "Vector.H"
#include "Tuple.H"

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
                int reflectDirection,
                int reflectCoordinate);

#include "NamespaceFooter.H"

#endif


