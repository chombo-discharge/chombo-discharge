#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _BINFABFACTORYI_H_
#define _BINFABFACTORYI_H_
#include "NamespaceHeader.H"

template <class T>
BinFabFactory<T>::~BinFabFactory()
{
}

template <class T>
void
BinFabFactory<T>::define(const RealVect&      a_meshSpacing,
                         const RealVect&      a_origin)
{
  m_mesh_spacing = a_meshSpacing;
  m_origin       = a_origin;
}

template <class T>
BinFab<T>*
BinFabFactory<T>::create(const Box &       a_box,
                         int               a_numComps,
                         const DataIndex & a_dit) const
{
  // [NOTE: a_ncomps and a_dit are ignored.]
  return new BinFab<T>(a_box, m_mesh_spacing, m_origin);
}

#include "NamespaceFooter.H"
#endif
