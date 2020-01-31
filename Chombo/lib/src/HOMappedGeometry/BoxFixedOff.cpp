#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxFixedOff.H"

#include "NamespaceHeader.H"

BoxFixedOff::BoxFixedOff(const Interval& a_fixedDims)
{
  m_fixedDims = a_fixedDims;
}

BoxFixedOff::~BoxFixedOff()
{
}

Box
BoxFixedOff::operator() (const Box& a_inputBox)
{
  Box returnBox(a_inputBox);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_fixedDims.contains(idir))
        { // Subtract low end in the fixed dimensions.
          returnBox.shift(idir, -a_inputBox.smallEnd(idir));
        }
      else
        { // Set to zero in the interpolated dimensions.
          returnBox.setRange(idir, 0); // default length 1 cell
        }
    }
  return returnBox;
}

#include "NamespaceFooter.H"
