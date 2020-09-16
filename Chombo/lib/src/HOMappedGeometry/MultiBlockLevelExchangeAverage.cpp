#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockLevelExchangeAverage.H"
#include "FourthOrderUtil.H"
#include "MBStencilIterator.H"
#include "MBVectorStencilIterator.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockLevelExchangeAverage::MultiBlockLevelExchangeAverage(const MultiBlockLevelGeom*  a_geomPtr,
                                                             int                         a_ghosts,
                                                             int                         a_order)
{
  define(a_geomPtr, a_ghosts, a_order);
}

/// destructor
MultiBlockLevelExchangeAverage::~MultiBlockLevelExchangeAverage()
{
  undefine();
}

void MultiBlockLevelExchangeAverage::define(const MultiBlockLevelGeom*  a_geomPtr,
                                           int                         a_ghosts,
                                           int                         a_order)
{
  undefine();
  m_type = IndexType::TheCellType();
  MultiBlockLevelExchange::define(a_geomPtr, a_ghosts, a_order);
  m_isDefined = true;
}

#include "NamespaceFooter.H"
