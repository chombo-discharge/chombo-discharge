/*!
  @file hyperboloid_one_if.cpp
  @brief Implementation of a one-sheeted hyperboloid surface
  @date Nov. 2017
  @author Robert Marskar
*/
  
#include "hyperboloid_one_if.H"

#include "CD_NamespaceHeader.H"

hyperboloid_one_if::hyperboloid_one_if(const RealVect& a_radii,
				       const RealVect& a_center,
				       const bool&     a_inside){
  m_radii  = a_radii;
  m_center = a_center;
  m_inside = a_inside;
  m_radii2 = m_radii*m_radii;
  m_sign   = RealVect::Unit;
  m_sign[SpaceDim-1] = -1.;
}

hyperboloid_one_if::hyperboloid_one_if(const hyperboloid_one_if& a_inputIF){
  m_radii  = a_inputIF.m_radii;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;
  m_radii2 = a_inputIF.m_radii2;
  m_sign   = a_inputIF.m_sign;
}

Real hyperboloid_one_if::value(const RealVect& a_point) const{
  Real retval;
  Real sum;

  // Compute the equation for the hyperboloid
  sum = 0.;
  for (int dir = 0; dir < SpaceDim; dir++){
    Real cur;
    cur = a_point[dir] - m_center[dir];
    sum += m_sign[dir]*cur*cur/m_radii2[dir];
  }

  // Sum is equal to -1 on the surface
  retval = sum - 1.0;

  // Change sign to change inside to outside
  if(!m_inside){
    retval = -retval;
  }

  return retval;
}

BaseIF* hyperboloid_one_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new hyperboloid_one_if(m_radii,m_center,m_inside));
}
#include "CD_NamespaceFooter.H"
