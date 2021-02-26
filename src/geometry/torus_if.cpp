/*!
  @file   torus_if.cpp
  @brief  Implementation of torus_if.H
  @date   Feb. 2021
  @author Robert Marskar
*/

#include "torus_if.H"


torus_if::torus_if(const RealVect a_center, const Real a_majorRadius, const Real a_minorRadius, const bool a_fluidInside){
  m_center      = a_center;
  m_majorRadius = a_majorRadius;
  m_minorRadius = a_minorRadius;
  m_fluidInside = a_fluidInside;
}


torus_if::torus_if(const torus_if& a_inputIF){
  m_center      = a_inputIF.m_center;
  m_majorRadius = a_inputIF.m_majorRadius;
  m_minorRadius = a_inputIF.m_minorRadius;
  m_fluidInside = a_inputIF.m_fluidInside;
}


torus_if::~torus_if(){

}

Real torus_if::value(const RealVect& a_point) const{
  const RealVect p  = a_point - m_center;

  Real radius = 0.0;
  for (int dir = 0; dir < 2; dir++){
    const Real cur = p[dir];
    radius += cur*cur;
  }
  radius = sqrt(radius) - m_majorRadius;
      

  Real retval = radius*radius;
#if CH_SPACEDIM==3
  retval += p[2]*p[2];
#endif

  retval = sqrt(retval) - m_minorRadius; // Positive outside. 

  if(!m_fluidInside){ // Make sure negative outside, if fluid is outside. 
    retval = -retval;
  }

  return retval;
}

BaseIF* torus_if::newImplicitFunction() const{
  return (BaseIF*) new torus_if(m_center, m_majorRadius, m_minorRadius, m_fluidInside);
}


