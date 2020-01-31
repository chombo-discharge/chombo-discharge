#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXFIXEDOFF_H_
#define _BOXFIXEDOFF_H_

#include "BoxLayout.H" // contains BaseTransform
#include "Interval.H"

#include "NamespaceHeader.H"

/// This class uncollapses the collapsed dimensions of a Box, and leaves other dimensions set to 0.

class BoxFixedOff : public BaseTransform
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  ///
  /** Constructor.
   */
  BoxFixedOff(const Interval& a_fixedDims);

  /** destructor
   */
  virtual ~BoxFixedOff();

  /** Required for BaseTransform.
      Returns the collapsed Box:
      the dimensions in a_fixedDims are set to a_fixedPt.
  */
  virtual Box operator() (const Box& a_inputBox);

protected:

  Interval m_fixedDims;
};

#include "NamespaceFooter.H"
#endif
