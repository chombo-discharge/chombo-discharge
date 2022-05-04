/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SignedDistanceDCELImplem.H
  @brief  Implementation of CD_SignedDistanceDCEL.H
  @author Robert Marskar
*/

#ifndef CD_SignedDistanceDCELImplem_H
#define CD_SignedDistanceDCELImplem_H

// Our includes
#include <CD_SignedDistanceDCEL.H>
#include <CD_NamespaceHeader.H>

template <class T>
SignedDistanceDCEL<T>::SignedDistanceDCEL(const std::shared_ptr<Mesh>& a_mesh,
                                          const bool                   a_flipInside,
                                          const Real                   a_zCoord)
{
  m_mesh       = a_mesh;
  m_zCoord     = a_zCoord;
  m_flipInside = a_flipInside;
}

template <class T>
SignedDistanceDCEL<T>::SignedDistanceDCEL(const SignedDistanceDCEL& a_object)
{
  m_mesh       = a_object.m_mesh;
  m_zCoord     = a_object.m_zCoord;
  m_flipInside = a_object.m_flipInside;
}

template <class T>
SignedDistanceDCEL<T>::~SignedDistanceDCEL()
{}

template <class T>
Real
SignedDistanceDCEL<T>::value(const RealVect& a_point) const
{

#if CH_SPACEDIM == 2
  EBGeometry::Vec3T<T> p(a_point[0], a_point[1], m_zCoord);
#else
  EBGeometry::Vec3T<T> p(a_point[0], a_point[1], a_point[2]);
#endif

  T retval = m_mesh->signedDistance(
    p); // Note that Dcel::mesh can return either positive or negative for outside, depending on the orientation of the input normals.

  if (m_flipInside) {
    retval = -retval;
  }

  return Real(retval);
}

template <class T>
BaseIF*
SignedDistanceDCEL<T>::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new SignedDistanceDCEL(*this));
}

#include <CD_NamespaceFooter.H>

#endif
