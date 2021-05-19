/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_HyperboloidTwoIF.cpp
  @brief  Implementation of CD_HyperboloidTwoIF.H
  @author Robert Marskar
*/

// Our includes
#include <CD_HyperboloidTwoIF.H>
#include <CD_NamespaceHeader.H>

HyperboloidTwoIF::HyperboloidTwoIF(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside){
  m_radii  = a_radii;
  m_center = a_center;
  m_inside = a_inside;
  m_radii2 = m_radii*m_radii;
  m_sign   = RealVect::Unit;
  m_sign[SpaceDim-1] = -1.;
}

HyperboloidTwoIF::HyperboloidTwoIF(const HyperboloidTwoIF& a_inputIF){
  m_radii  = a_inputIF.m_radii;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;
  m_radii2 = a_inputIF.m_radii2;
  m_sign   = a_inputIF.m_sign;
}

Real HyperboloidTwoIF::value(const RealVect& a_point) const{
  Real retval;
  Real sum;

  for (int dir = 0; dir < SpaceDim; dir++){
    Real cur;
    cur = a_point[dir] - m_center[dir];
    sum += m_sign[dir]*cur*cur/m_radii2[dir];
  }

  retval = sum + 1.0;

  if(!m_inside){
    retval = -retval;
  }

  return retval;
}

BaseIF* HyperboloidTwoIF::newImplicitFunction() const{
  return static_cast<BaseIF*> (new HyperboloidTwoIF(m_radii,m_center,m_inside));
}

#include <CD_NamespaceFooter.H>
