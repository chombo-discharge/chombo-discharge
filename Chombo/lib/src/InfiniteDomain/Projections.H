#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Projections.H
// petermc, 2 Sep 2004

#ifndef _PROJECTIONS_H_
#define _PROJECTIONS_H_

#include "FArrayBox.H"

#include "NamespaceHeader.H"

///
/**
   Returns the "outer" coarsening of a NODE-centered box:
   coarsenOuter(((l1, l2, l3):(u1, u2, u3)), r) ==
   (floor(l1/r, l2/r, l3/r):ceil(u1/r, u2/r, u3/r))

   If r divides all of l1, l2, l3, u1, u2, u3 then we have simply
   coarsenOuter(((l1, l2, l3):(u1, u2, u3)), r) ==
   coarsen(((l1, l2, l3):(u1, u2, u3)), r)
   == ((l1/r, l2/r, l3/r):(u1/r, u2/r, u3/r))
 */
Box coarsenOuter(/// NODE-centered fine box
                 const Box&   a_bx,
                 /// refinement ratio
                 int          a_refRatio);

///
/**
   Returns the "inner" coarsening of a NODE-centered box:
   coarsenInner(((l1, l2, l3):(u1, u2, u3)), r) ==
   (ceil(l1/r, l2/r, l3/r):floor(u1/r, u2/r, u3/r))

   If r divides all of l1, l2, l3, u1, u2, u3 then we have simply
   coarsenInner(((l1, l2, l3):(u1, u2, u3)), r) ==
   coarsen(((l1, l2, l3):(u1, u2, u3)), r)
   == ((l1/r, l2/r, l3/r):(u1/r, u2/r, u3/r))
 */
Box coarsenInner(/// NODE-centered fine box
                 const Box&   a_bx,
                 /// refinement ratio
                 int          a_refRatio);

///
/**
   Projects a fine FAB to a coarse FAB by
   a_coarse(ivCoarse) = a_fine(ivCoarse * a_refRatio + a_shift)
   where a_shift is optional, defaulting to zero.

   Data are projected from nodes in the intersection of the box of a_coarse
   and the inner coarsening (see coarsenInner above) by a_refRatio of
   the box of a_fine shifted by -a_shift.
   Both FABs must have only one component.
 */
void projectToCoarse(/// coarse FAB result of projection
                     FArrayBox&         a_coarse,
                     /// fine FAB of one component to be projected
                     const FArrayBox&   a_fine,
                     /// refinement ratio
                     int                a_refRatio,
                     /// amount to shift fine FAB before projecting
                     const IntVect&     a_shift = IntVect::Zero);

///
/**
   Projects a coarse FAB to a fine FAB by
   a_coarse(ivCoarse * a_refRatio + a_shift) = a_coarse(ivCoarse)
   where a_shift is optional, defaulting to zero.

   Data are projected from nodes in the intersection of the box of a_coarse
   and the inner coarsening (see coarsenInner above) by a_refRatio of
   the box of a_fine shifted by -a_shift.
   Both FABs must have only one component.
 */
void projectToFine(/// fine FAB result of projection
                   FArrayBox&         a_fine,
                   /// coarse FAB of one component to be projected
                   const FArrayBox&   a_coarse,
                   /// refinement ratio
                   int                a_refRatio,
                   /// amount to shift fine FAB after projecting
                   const IntVect&     a_shift = IntVect::Zero);

#include "NamespaceFooter.H"

#endif
