#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXCOLLAPSER_H_
#define _BOXCOLLAPSER_H_

#include "BoxLayout.H" // contains BaseTransform
#include "Interval.H"

#include "NamespaceHeader.H"

/// This class collapses the "fixed" dimensions of a Box to its small end.

class BoxCollapser : public BaseTransform
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  ///
  /** Constructor.
   */
  BoxCollapser(const Interval& a_fixedDims);

  /** destructor
   */
  virtual ~BoxCollapser();

  /** Required for BaseTransform.
      Returns the collapsed Box:
      the range in a_fixedDims is set to a_inputBox.smallEnd() only.
  */
  virtual Box operator() (const Box& a_inputBox);

protected:

  Interval m_fixedDims;
};

#include "NamespaceFooter.H"
#endif
