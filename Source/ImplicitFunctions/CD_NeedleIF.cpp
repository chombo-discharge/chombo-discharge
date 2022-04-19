/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.cpp
  @brief  Implementation of CD_NeedleIF.H
  @author Fanny Skirbekk
*/

// Our includes
#include <CD_NeedleIF.H>
#include <CD_CylinderSdf.H>
#include <CD_NamespaceHeader.H>
#include <EBGeometry_AnalyticDistanceFunctions.hpp> // can I just do like this?

NeedleIF::NeedleIF(const RealVect& a_centerTipSide, const Realvect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipRadius, const Real& a_angle){

  constexpr Real pi = 3.14159265358979323846;
  const Real tipLength = a_radius/std::tan(angle*pi/180);

  const RealVect axis    = (a_centerTipSide - a_centerBack);
  const RealVect axisVec = axis/axis.vectorLength(); 

  // find new center for where the cylinder and cone should meet.
  const RealVect c = a_centerTipSide - axisVec*tipLength;

  // Build the needle-parts
  Vector<BaseIF*> isects;
  isects.push_back(static_cast<BaseIF*> (new CylinderSdf(c, a_centerBack, a_radius, a_fluidInside)));
  isects.push_back(static_cast<BaseIF*> (new ConeSDF(a_centerTipSide, tipLength, a_angle, false));

  // Build the needle
  m_baseif = RefCountedPtr<BaseIF>(new IntersectionIF(isects));

  // Delete everything we have allocated so far
  for(int i = 0; i < isects.size(); ++i{
    delete isects[i];
  }
}

NeedleIF::NeedleIF(const NeedleIF& a_inputIF){
  this->m_baseif = a_inputIF.m_baseif;
}

Real NeedleIF::value(const RealVect& a_point) const{
  return m_baseif->value(a_point);
}

BaseIF* NeedleIF::newImplicitFunction() const{
  return static_cast<BaseIF*> (new NeedleIF(*this));
}

#include <CD_NamespaceFooter.H>
