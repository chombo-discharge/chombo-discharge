#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GROWINBLOCK_H_
#define _GROWINBLOCK_H_

// class BaseTransform is declared in BoxLayout.H
#include "BoxLayout.H"
#include "LevelGridMetrics.H"

#include "NamespaceHeader.H"

/// Class that is given a base box in a MultiBlockCoordSys and returns the chunk of its halo to be filled in from a specified source block with a specified offset.
class GrowInBlock : public BaseTransform
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  ///
  /** Constructor.
   */
  GrowInBlock(
              LevelGridMetrics* a_levelGridMetricsPtr,
              int a_ghost);

  /** destructor
   */
  virtual ~GrowInBlock();

  /** Required for BaseTransform.
      Returns the Box consisting of the input box grown by m_radius,
      but cut off at physical boundaries.
  */
  virtual Box operator() (const Box& a_inputBox);

protected:

  /// width of ghost layer
  int m_ghost;

  LevelGridMetrics* m_levelGridMetricsPtr;
};

#include "NamespaceFooter.H"

#endif // include guard
