/*!
  @file perlin_sphere_if.H
  @brief Implementation of perlin_sphere_if.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "perlin_sphere_if.H"

#include "CD_NamespaceHeader.H"
  
perlin_sphere_if::perlin_sphere_if(const Real&     a_rad,
				   const RealVect& a_center,
				   const bool&     a_inside,
				   const Real&     a_noiseAmp,
				   const RealVect& a_noiseFreq,
				   const Real&     a_persistence,
				   const int&      a_octaves,
				   const bool&     a_reseed){

  //
  m_rad    = a_rad;
  m_center = a_center;
  m_inside = a_inside;
  m_perlinIF = RefCountedPtr<BaseIF> (new perlin_if(a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed));
}

//
perlin_sphere_if::perlin_sphere_if(const perlin_sphere_if& a_inputIF){

  //
  m_rad      = a_inputIF.m_rad;
  m_center   = a_inputIF.m_center;
  m_inside   = a_inputIF.m_inside;
  m_perlinIF = a_inputIF.m_perlinIF;
}

//
perlin_sphere_if::~perlin_sphere_if(){
}

//
Real perlin_sphere_if::value(const RealVect& a_pos) const {

  RealVect v;
  Real retval;

  const RealVect pos = a_pos - m_center;
  
  // Get noise on the circle/sphere
#if CH_SPACEDIM == 2
  const Real theta = atan2(pos[0], pos[1]);
  const Real x     = m_rad*sin(theta);
  const Real y     = m_rad*cos(theta);
  
  v = RealVect(x,y);
#elif CH_SPACEDIM == 3
  const Real xy    = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  const Real theta = atan2(xy, pos[2]);
  const Real phi   = atan2(pos[1], pos[0]);
  const Real x     = m_rad*sin(theta)*sin(phi);
  const Real y     = m_rad*sin(theta)*cos(phi);
  const Real z     = m_rad*cos(theta);
  
  v = RealVect(x,y,z);
#endif

  // Random radius
  const Real R = m_rad + m_perlinIF->value(v);

  // Value function
  const Real dist2 = pos.vectorLength()*pos.vectorLength() - R*R;
  if(dist2 > 0.){
    retval = sqrt(dist2);
  }
  else{
    retval = -sqrt(-dist2);
  }

  // Switch inside to outside
  if(!m_inside){
    retval = -retval;
  }
  
  return retval;
}

//
BaseIF* perlin_sphere_if::newImplicitFunction() const {
  return static_cast<BaseIF*> (new perlin_sphere_if(*this));
}
#include "CD_NamespaceFooter.H"
