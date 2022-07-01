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

// CHombo includes
#include <UnionIF.H>
#include <SmoothIntersection.H>

// Our includes
#include <CD_NeedleIF.H>
#include <CD_CylinderSdf.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using Vec3 = EBGeometry::Vec3T<Real>;

NeedleIF::NeedleIF(const RealVect& a_centerTipSide, const RealVect& a_centerBack, const Real& a_radius, const bool& a_fluidInside, const Real& a_tipRadius, const double& a_angle){
  m_tipRadius = {0, a_tipRadius, 0}; 
  
  constexpr Real pi = 3.14159265358979323846;
  const double tipLength = a_radius/std::tan(a_angle*pi/180);

  const RealVect axis    = (a_centerTipSide - a_centerBack);
  const RealVect axisVec = axis/axis.vectorLength(); 

  // find new center for where the cylinder and cone should meet.
  const RealVect c = a_centerTipSide - tipLength*axisVec;

#if CH_SPACEDIM==2
    const Vec3 centerT(a_centerTipSide[0], a_centerTipSide[1], 0.0);
#else
    const Vec3 centerT(a_centerTipSide[0], a_centerTipSide[1], a_centerTipSide[2]);
#endif

  // Build the needle-parts
  Vector<BaseIF*> isects;
  isects.push_back(static_cast<BaseIF*> (new CylinderSdf(c, a_centerBack, a_radius/2, a_fluidInside)));
  //flipinside=true for cone since EBGeometry and Chombo has opposing sign conventions/logic regarding the flipinside..
  auto cone = std::make_shared<EBGeometry::ConeSDF<Real>>(centerT, tipLength, a_angle,true);
  //for debugging..:
  /*
  std::cout << "axisVec:" << axisVec[0] << "," << axisVec[1] << "," << axisVec[2] << std::endl;
  std::cout << "centerT:" << centerT[0] << "," << centerT[1] << "," << centerT[2] << std::endl;
  std::cout << "c:" << c[0] << "," << c[1] << "," << c[2] << std::endl;
  std::cout << "tiplength:" << tipLength << std::endl;
  */
  cone->rotate(90,0);
  isects.push_back(static_cast<BaseIF*> (new EBGeometryIF(cone, false)));

  // Build the needle
  m_baseif = RefCountedPtr<BaseIF>(new SmoothIntersection(isects, 0.01));

  // Delete everything we have allocated so far
  for(int i = 0; i < isects.size(); ++i){
    delete isects[i];
  }
}

NeedleIF::NeedleIF(const NeedleIF& a_inputIF){
  this->m_baseif = a_inputIF.m_baseif;
}

Real NeedleIF::value(const RealVect& a_point) const{
  // tried to "round" of the point
  return m_baseif->value(a_point - m_tipRadius);
}

BaseIF* NeedleIF::newImplicitFunction() const{
  return static_cast<BaseIF*> (new NeedleIF(*this));
}

#include <CD_NamespaceFooter.H>
