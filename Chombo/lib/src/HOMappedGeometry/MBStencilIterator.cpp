#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBStencilIterator.H"
#include "NamespaceHeader.H"

void MBStencilIterator::define(const MBStencil& a_stencil)
{
  m_stencil = a_stencil;
  m_indCurrent = 0;
  m_indLo = 0;
  m_indHi = a_stencil.size() - 1;
}

void MBStencilIterator::setStencil(const MBStencil& a_stencil)
{
  define(a_stencil);
}

#include "NamespaceFooter.H"
