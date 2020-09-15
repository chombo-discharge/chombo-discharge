#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxCollapser.H"

#include "NamespaceHeader.H"

BoxCollapser::BoxCollapser(const Interval& a_fixedDims)
{
  m_fixedDims = a_fixedDims;
}

BoxCollapser::~BoxCollapser()
{
}

Box
BoxCollapser::operator() (const Box& a_inputBox)
{
  Box returnBox(a_inputBox);
  for (int idir = m_fixedDims.begin(); idir <= m_fixedDims.end(); idir++)
    { // Collapse range in direction dir a_inputsBox.smallEnd(dir).
      int lo = a_inputBox.smallEnd(idir);
      returnBox.setRange(idir, lo); // default length 1 cell
    }
  return returnBox;
}

#include "NamespaceFooter.H"
