/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGeometryIF.cpp
  @brief  Implementation of CD_EBGeometryIF.H
  @author Robert marskar
*/

#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

EBGeometryIF::EBGeometryIF() {
  this->m_sdf        = nullptr;
  this->m_flipInside = false;
}

EBGeometryIF::EBGeometryIF(const std::shared_ptr<EBGeometry::SignedDistanceFunction<Real> >& a_sdf, const bool a_flipInside) {
  this->m_sdf        = a_sdf;
  this->m_flipInside = a_flipInside;
}

EBGeometryIF::EBGeometryIF(const EBGeometryIF& a_inputIF) {
  this->m_sdf        = a_inputIF.m_sdf;
  this->m_flipInside = a_inputIF.m_flipInside;  
}

EBGeometryIF::~EBGeometryIF() {
  m_sdf = nullptr;
}

Real EBGeometryIF::value(const RealVect& a_point) const {
#if CH_SPACEDIM==2
  EBGeometry::Vec3T<Real> p(a_point[0], a_point[1], 0.0);
#else
  EBGeometry::Vec3T<Real> p(a_point[0], a_point[1], a_point[2]);
#endif

  Real ret = Real(m_sdf->signedDistance(p));

  if(m_flipInside) {
    ret = -ret;
  }

  return ret;
}

BaseIF* EBGeometryIF::newImplicitFunction() const {
  return (BaseIF*) (new EBGeometryIF(*this));
}

#include <CD_NamespaceFooter.H>
