#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBStencilElement.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
MBStencilElement::MBStencilElement()
{
  m_block = -1;
}

// ---------------------------------------------------------
MBStencilElement::MBStencilElement(const IntVect&  a_cell,
                                   int             a_block,
                                   Real            a_weight)
{
  define(a_cell, a_block, a_weight);
}


// ---------------------------------------------------------
void
MBStencilElement::define(const IntVect&  a_cell,
                         int             a_block,
                         Real            a_weight)
{
  m_cell = a_cell;
  m_block = a_block;
  m_weight = a_weight;
}

// ---------------------------------------------------------
bool
MBStencilElement::operator==(const MBStencilElement& a_input) const
{
  bool retval = ( (a_input.m_cell == m_cell) &&
                  (a_input.m_block == m_block) &&
                  (a_input.m_weight == m_weight) );
  return retval;
}

#include "NamespaceFooter.H"
