#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELCOPIERAVERAGE_H_
#define _MULTIBLOCKLEVELCOPIERAVERAGE_H_

#include "MultiBlockLevelCopier.H"
#include "VectorAvgTransformation.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelCopierAverage stores topological and geometric information
/**
   MultiBlockLevelCopierAverage stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys,
   for copying to another layout, on cell-averaged data.
*/
class MultiBlockLevelCopierAverage : public MultiBlockLevelCopier
{


public:

  /// default constructor
  MultiBlockLevelCopierAverage() : MultiBlockLevelCopier()
  {
  }

  /// destructor
  virtual ~MultiBlockLevelCopierAverage();

  /// full constructor
  MultiBlockLevelCopierAverage(const MultiBlockLevelGeom*  a_geomPtr,
                               const BoxLayout&            a_dstLayout,
                               int                         a_order);

  MultiBlockLevelCopierAverage(const MultiBlockLevelGeom*  a_geomPtr,
                               const DisjointBoxLayout&    a_dstDisjointLayout,
                               const IntVect&              a_ghostVect,
                               int                         a_order);

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      const BoxLayout&            a_dstLayout,
                      int                         a_order);

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      const DisjointBoxLayout&    a_dstDisjointLayout,
                      const IntVect&              a_ghostVect,
                      int                         a_order);

protected:

};

#include "NamespaceFooter.H"
#endif
