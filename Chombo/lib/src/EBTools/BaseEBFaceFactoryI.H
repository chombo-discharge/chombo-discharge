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

#ifndef _BASEEBFACEFACTORYI_H_
#define _BASEEBFACEFACTORYI_H_
#include "NamespaceHeader.H"

/***************/
/***************/
template <class T>
BaseEBFaceFactory<T>::~BaseEBFaceFactory()
{
}
/***************/
/***************/
template <class T>
BaseEBFaceFactory<T>::BaseEBFaceFactory(const EBISLayout& a_ebisl, int a_idir)
{
  m_ebisl = a_ebisl;
  m_idir = a_idir;
}
/***************/
/***************/
template <class T>
BaseEBFaceFAB<T>*
BaseEBFaceFactory<T>::create(const Box& a_box, int a_ncomps,
                             const DataIndex& a_dit) const
{
  return new BaseEBFaceFAB<T>(m_ebisl[a_dit], a_box, m_idir, a_ncomps);
}
/***************/
/***************/

#include "NamespaceFooter.H"
#endif
