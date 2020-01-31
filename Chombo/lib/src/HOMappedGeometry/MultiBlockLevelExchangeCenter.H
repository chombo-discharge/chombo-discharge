#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELEXCHANGECENTER_H_
#define _MULTIBLOCKLEVELEXCHANGECENTER_H_

#include "MultiBlockLevelExchange.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelExchangeCenter stores topological and geometric information for cell-centered data
/**
   MultiBlockLevelExchangeCenter stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys,
   on cell-centered data.
*/
class MultiBlockLevelExchangeCenter : public MultiBlockLevelExchange
{

public:

  /// default constructor
  MultiBlockLevelExchangeCenter() : MultiBlockLevelExchange()
  {
  }

  /// destructor
  virtual ~MultiBlockLevelExchangeCenter();

  /// full constructor
  MultiBlockLevelExchangeCenter(const MultiBlockLevelGeom*  a_geomPtr,
                                int                         a_ghosts,
                                int                         a_order);

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      int                         a_ghosts,
                      int                         a_order);

protected:

};

#include "NamespaceFooter.H"
#endif
