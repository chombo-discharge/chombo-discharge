/*!
  @file hyperboloid_two_if.cpp
  @brief Implementation of a one-sheeted hyperboloid surface
  @date Nov. 2015
  @author Robert Marskar
*/
  
#include "hyperboloid_two_if.H"

hyperboloid_two_if::hyperboloid_two_if(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside){
  m_radii  = a_radii;
  m_center = a_center;
  m_inside = a_inside;
  m_radii2 = m_radii*m_radii;
  m_sign   = RealVect::Unit;
  m_sign[SpaceDim-1] = -1.;
}

hyperboloid_two_if::hyperboloid_two_if(const hyperboloid_two_if& a_inputIF){
  m_radii  = a_inputIF.m_radii;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;
  m_radii2 = a_inputIF.m_radii2;
  m_sign   = a_inputIF.m_sign;
}

Real hyperboloid_two_if::value(const RealVect& a_point) const{
  Real retval;
  Real sum;

  //
  for (int dir = 0; dir < SpaceDim; dir++){
    Real cur;
    cur = a_point[dir] - m_center[dir];
    sum += m_sign[dir]*cur*cur/m_radii2[dir];
  }

  //
  retval = sum + 1.0;

  //
  if(!m_inside){
    retval = -retval;
  }

  return retval;
}

BaseIF* hyperboloid_two_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new hyperboloid_two_if(m_radii,m_center,m_inside));
}
