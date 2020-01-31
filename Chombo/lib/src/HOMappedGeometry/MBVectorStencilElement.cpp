#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBVectorStencilElement.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
MBVectorStencilElement::MBVectorStencilElement()
{
  m_block = -1;
}

// ---------------------------------------------------------
MBVectorStencilElement::MBVectorStencilElement(const IntVect&  a_cell,
                                               int             a_block,
                                               Tuple<Real, SpaceDim*SpaceDim> a_weight)
{
  define(a_cell, a_block, a_weight);
}


// ---------------------------------------------------------
MBVectorStencilElement::MBVectorStencilElement(const IntVect&  a_cell,
                                               int             a_block)
{
  define(a_cell, a_block);
}


// ---------------------------------------------------------
void
MBVectorStencilElement::define(const IntVect&  a_cell,
                               int             a_block,
                               Tuple<Real, SpaceDim*SpaceDim> a_weight)
{
  m_cell = a_cell;
  m_block = a_block;
  m_weight = a_weight;
}

// ---------------------------------------------------------
void
MBVectorStencilElement::define(const IntVect&  a_cell,
                               int             a_block)
{
  m_cell = a_cell;
  m_block = a_block;
  for (int i = 0; i < SpaceDim*SpaceDim; i++)
    {
      m_weight[i] = 0.;
    }
}

// ---------------------------------------------------------
bool
MBVectorStencilElement::operator==(const MBVectorStencilElement& a_input) const
{
  bool sameCell = (a_input.m_cell == m_cell);
  bool sameBlock = (a_input.m_block == m_block);
  bool sameWeight = true;
  for (int i = 0; i < SpaceDim*SpaceDim; i++)
    {
      sameWeight &= (a_input.m_weight[i] == m_weight[i]);
    }
  bool sameAll = sameCell && sameBlock && sameWeight;

  return sameAll;
}

#include "NamespaceFooter.H"
