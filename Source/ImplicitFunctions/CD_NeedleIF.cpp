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
#include <UnionIF.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using Vec3 = EBGeometry::Vec3T<Real>;

NeedleIF::NeedleIF(const RealVect& a_centerTipSide, const RealVect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipRadius, const double& a_angle){

  constexpr Real pi = 3.14159265358979323846;
  const double tipLength = a_radius/std::tan(a_angle*pi/180);

  const RealVect axis    = (a_centerTipSide - a_centerBack);
  const RealVect axisVec = axis/axis.vectorLength(); 

  // find new center for where the cylinder and cone should meet.
  const RealVect c = a_centerTipSide - tipLength*axisVec;

#if CH_SPACEDIM==2
    const Vec3 centerT(0.0, 0.0, 0.0);
#else
    const Vec3 centerT(0.0, 0.0, 0.0);
#endif

  // Build the needle-parts
  Vector<BaseIF*> isects;
  //isects.push_back(static_cast<BaseIF*> (new CylinderSdf(c, a_centerBack, a_radius, a_fluidInside)));
  auto cone = std::make_shared<EBGeometry::ConeSDF<Real>>(centerT, 0.5, a_angle, false);
  std::cout << c << " " <<  tipLength << " " << a_angle << "\n"; 
  cone->rotate(90,0);
  isects.push_back(static_cast<BaseIF*> (new EBGeometryIF(cone, false)));

  // Build the needle
  m_baseif = RefCountedPtr<BaseIF>(new UnionIF(isects));

  // Delete everything we have allocated so far
  for(int i = 0; i < isects.size(); ++i){
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
