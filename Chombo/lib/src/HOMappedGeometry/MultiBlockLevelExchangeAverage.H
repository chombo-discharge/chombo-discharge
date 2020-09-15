#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELEXCHANGEAVERAGE_H_
#define _MULTIBLOCKLEVELEXCHANGEAVERAGE_H_

#include "MultiBlockLevelExchange.H"
#include "VectorAvgTransformation.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelExchangeAverage stores topological and geometric information for cell-averaged data
/**
   MultiBlockLevelExchangeAverage stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys,
   on cell-averaged data.
*/
class MultiBlockLevelExchangeAverage : public MultiBlockLevelExchange
{

public:

  /// default constructor
  MultiBlockLevelExchangeAverage() : MultiBlockLevelExchange()
  {
  }

  /// destructor
  virtual ~MultiBlockLevelExchangeAverage();

  /// full constructor
  MultiBlockLevelExchangeAverage(const MultiBlockLevelGeom*  a_geomPtr,
                                int                         a_ghosts,
                                int                         a_order);

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      int                         a_ghosts,
                      int                         a_order);

protected:

};

#include "NamespaceFooter.H"
#endif
