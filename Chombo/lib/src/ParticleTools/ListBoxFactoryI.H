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

#ifndef _LISTBOXFACTORYI_H_
#define _LISTBOXFACTORYI_H_
#include "ListBoxFactory.H"
#include "NamespaceHeader.H"

template <class T>
ListBoxFactory<T>::ListBoxFactory(){}

template <class T>
ListBoxFactory<T>::ListBoxFactory(const RealVect& a_meshSpacing,
                                  const RealVect& a_origin)
{
  define(a_meshSpacing, a_origin);
}

template <class T>
ListBoxFactory<T>::ListBoxFactory(const Real& a_meshSpacing)
{
  define(a_meshSpacing * RealVect::Unit, RealVect::Zero);
}

template <class T>
void
ListBoxFactory<T>::define(const RealVect& a_meshSpacing,
                          const RealVect& a_origin)
{
  m_mesh_spacing = a_meshSpacing;
  m_origin = a_origin;
}

template <class T>
ListBox<T>* 
ListBoxFactory<T>::create(const Box& a_box, int a_numComps,
			  const DataIndex& a_datInd) const
{
  return new ListBox<T>(a_box, m_mesh_spacing, m_origin);
}

/// Destructor
template <class T>
ListBoxFactory<T>::~ListBoxFactory()
{}

template <class T>
RealVect
ListBoxFactory<T>::meshSpacing() const
{
  return m_mesh_spacing;
}

template <class T>
RealVect
ListBoxFactory<T>::origin() const
{
  return m_origin;
}

#include "NamespaceFooter.H"
#endif   
