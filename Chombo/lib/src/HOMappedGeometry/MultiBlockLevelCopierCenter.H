#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELCOPIERCENTER_H_
#define _MULTIBLOCKLEVELCOPIERCENTER_H_

#include "MultiBlockLevelCopier.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelCopierCenter stores topological and geometric information
/**
   MultiBlockLevelCopierCenter stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys,
   for copying to another layout, on cell-centered data.
*/
class MultiBlockLevelCopierCenter : public MultiBlockLevelCopier
{


public:

  /// default constructor
  MultiBlockLevelCopierCenter() : MultiBlockLevelCopier()
  {
  }

  /// destructor
  virtual ~MultiBlockLevelCopierCenter();

  /// full constructor
  MultiBlockLevelCopierCenter(const MultiBlockLevelGeom*  a_geomPtr,
                              const BoxLayout&            a_dstLayout,
                              int                         a_order);

  MultiBlockLevelCopierCenter(const MultiBlockLevelGeom*  a_geomPtr,
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
