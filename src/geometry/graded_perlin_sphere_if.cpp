/*!
  @file   graded_perlin_sphere_if.cpp
  @brief  Implementation of graded_perlin_sphere_if.H
  @author Robert Marskar
  @date   Oct. 2017
*/

#include <math.h>

#include "graded_perlin_sphere_if.H"

#include "CD_NamespaceHeader.H"

graded_perlin_sphere_if::graded_perlin_sphere_if(const Real&     a_rad,
						 const RealVect& a_center,
						 const bool&     a_inside,
						 const Real&     a_noiseAmp,
						 const RealVect& a_noiseFreq,
						 const Real&     a_persistence,
						 const int&      a_octaves,
						 const bool&     a_reseed) : perlin_sphere_if(a_rad,
											      a_center,
											      a_inside,
											      a_noiseAmp,
											      a_noiseFreq,
											      a_persistence,
											      a_octaves,
											      a_reseed){
}

//
graded_perlin_sphere_if::graded_perlin_sphere_if(const graded_perlin_sphere_if& a_inputIF) : perlin_sphere_if(a_inputIF){
}

//
graded_perlin_sphere_if::~graded_perlin_sphere_if(){
}

//
Real graded_perlin_sphere_if::value(const RealVect& a_pos) const {

  //
  RealVect v;
  Real theta;
  
  const RealVect pos = a_pos - m_center;
  
  // Get noise on the circle/sphere
#if CH_SPACEDIM == 2
  theta        = atan2(pos[0], pos[1]);
  const Real x = m_rad*sin(theta);
  const Real y = m_rad*cos(theta);
  
  v = RealVect(x,y);
#elif CH_SPACEDIM == 3
  theta          = atan2(sqrt(pos[0]*pos[0] + pos[1]*pos[1]), pos[2]);
  const Real phi = atan2(pos[1], pos[0]);
  const Real x   = m_rad*sin(theta)*sin(phi);
  const Real y   = m_rad*sin(theta)*cos(phi);
  const Real z   = m_rad*cos(theta);
  
  v = RealVect(x,y,z);
#endif

  // Random radius using Perlin function
  const Real dTheta = abs(theta);
  const Real factor = (dTheta < M_PI/2.) ? 1 - dTheta*2./M_PI : 0;
  const Real R = m_rad + m_perlinIF->value(v)*factor;

  // Value function
  Real retval;
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
BaseIF* graded_perlin_sphere_if::newImplicitFunction() const {
  return static_cast<BaseIF*> (new graded_perlin_sphere_if(*this));
}
#include "CD_NamespaceFooter.H"
