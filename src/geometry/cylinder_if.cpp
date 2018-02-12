/*!
  @file cylinder_if.cpp
  @brief Implementation of cylinder_if.H
  @date Nov. 2017
  @author Robert Marskar
*/

#include <PolyGeom.H>

#include "cylinder_if.H"

cylinder_if::cylinder_if(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_inside){
  m_center1 = a_center1;
  m_center2 = a_center2;
  m_length  = sqrt(PolyGeom::dot(m_center2-m_center1,m_center2-m_center1));
  m_radius  = a_radius;
  m_inside  = a_inside;
}

cylinder_if::cylinder_if(const cylinder_if& a_inputIF){
  m_center1 = a_inputIF.m_center1;
  m_center2 = a_inputIF.m_center2;
  m_length  = a_inputIF.m_length;
  m_radius  = a_inputIF.m_radius;
  m_inside  = a_inputIF.m_inside;
}

Real cylinder_if::value(const RealVect& a_point) const{
  
  const RealVect cylTop  = m_center2 - m_center1;
  const RealVect cylAxis = cylTop/sqrt(PolyGeom::dot(cylTop,cylTop));

  // Translate cylinder center to origin and do all calculations for "positive half" of the cylinder
  const RealVect newPoint = a_point - m_center1 - (m_center2-m_center1)*0.5;
  const Real paraComp     = PolyGeom::dot(newPoint,cylAxis);
  const RealVect orthoVec = newPoint - paraComp*cylAxis;
  const Real orthoComp    = sqrt(PolyGeom::dot(orthoVec,orthoVec));

  const Real f = orthoComp - m_radius;
  const Real g = abs(paraComp)  - 0.5*m_length;

  Real retval = -1;
  if(f <= 0. && g <= 0.){      // Point lies within the cylinder. Either short end or wall is closest point. 
    retval = f <= g ? f : g;
  }
  else if(f <= 0. && g >= 0.){ // Point lies inside radius but outside length. Short end is the closest point
    retval = g;
  }
  else if(f > 0. && g <= 0.){  // Point lies outside radius but inside length. Cylinder wall is the closest point. 
    retval = f;
  }
  else if(f > 0. && g >  0.){  // Point lies outside both. The cylinder corner is the closest point.
    retval = sqrt(f*f + g*g);
  }

  // Change sign if called for
  if(!m_inside){
    retval = -retval;
  }

  return retval;
}

BaseIF* cylinder_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new cylinder_if(m_center1,m_center2,m_radius,m_inside));
}

