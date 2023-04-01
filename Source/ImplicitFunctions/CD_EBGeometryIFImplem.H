/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGeometryIFImplem.H
  @brief  Implementation of CD_EBGeometryIF.H
  @author Robert marskar
*/

#ifndef CD_EBGeometryIFImplem_H
#define CD_EBGeometryIFImplem_H

#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

template <typename T>
EBGeometryIF<T>::EBGeometryIF()
{
  this->m_sdf        = nullptr;
  this->m_flipInside = false;
  this->m_zCoord     = 0.0;
}

template <typename T>
EBGeometryIF<T>::EBGeometryIF(const std::shared_ptr<EBGeometry::ImplicitFunction<T>>& a_sdf,
                              const bool                                              a_flipInside,
                              const Real                                              a_zCoord)

{
  this->m_sdf        = a_sdf;
  this->m_flipInside = a_flipInside;
  this->m_zCoord     = a_zCoord;
}

template <typename T>
EBGeometryIF<T>::EBGeometryIF(const EBGeometryIF& a_inputIF)
{
  this->m_sdf        = a_inputIF.m_sdf;
  this->m_flipInside = a_inputIF.m_flipInside;
  this->m_zCoord     = a_inputIF.m_zCoord;
}

template <typename T>
EBGeometryIF<T>::~EBGeometryIF()
{
  m_sdf = nullptr;
}

template <typename T>
Real
EBGeometryIF<T>::value(const RealVect& a_point) const
{
#if CH_SPACEDIM == 2
  EBGeometry::Vec3T<T> p(a_point[0], a_point[1], m_zCoord);
#else
  EBGeometry::Vec3T<T> p(a_point[0], a_point[1], a_point[2]);
#endif

  Real ret = Real(m_sdf->value(p));

  if (m_flipInside) {
    ret = -ret;
  }

  return ret;
}

template <typename T>
BaseIF*
EBGeometryIF<T>::newImplicitFunction() const
{
  return (BaseIF*)(new EBGeometryIF(*this));
}

#include <CD_NamespaceFooter.H>

#endif
