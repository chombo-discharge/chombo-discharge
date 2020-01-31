#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockLevelCopierAverage.H"
#include "FourthOrderUtil.H"
#include "CH_Timer.H"
#include "MBStencilIterator.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockLevelCopierAverage::MultiBlockLevelCopierAverage(const MultiBlockLevelGeom*  a_geomPtr,
                                                         const BoxLayout&            a_dstLayout,
                                                         int                         a_order)
{
  define(a_geomPtr, a_dstLayout, a_order);
}


MultiBlockLevelCopierAverage::MultiBlockLevelCopierAverage(const MultiBlockLevelGeom*  a_geomPtr,
                                                           const DisjointBoxLayout&    a_dstDisjointLayout,
                                                           const IntVect&              a_ghostVect,
                                                           int                         a_order)
{
  define(a_geomPtr, a_dstDisjointLayout, a_ghostVect, a_order);
}


/// destructor
MultiBlockLevelCopierAverage::~MultiBlockLevelCopierAverage()
{
  undefine();
}


void
MultiBlockLevelCopierAverage::define(const MultiBlockLevelGeom*  a_geomPtr,
                                     const DisjointBoxLayout&    a_dstDisjointLayout,
                                     const IntVect&              a_ghostVect,
                                     int                         a_order)
{
  undefine();
  BoxLayout dstLayout;
  dstLayout.deepCopy(a_dstDisjointLayout);
  dstLayout.close();
  dstLayout.grow(a_ghostVect);
  define(a_geomPtr, dstLayout, a_order);
}


void
MultiBlockLevelCopierAverage::define(const MultiBlockLevelGeom*  a_geomPtr,
                                     const BoxLayout&            a_dstLayout,
                                     int                         a_order)
{
  undefine();
  m_type = IndexType::TheCellType();
  MultiBlockLevelCopier::define(a_geomPtr, a_dstLayout, a_order);
  m_isDefined = true;
}

#include "NamespaceFooter.H"
