/*!
  @file    profile_cylinder_if.cpp
  @brief   Implementation of profile_cylinder_if.H
  @author  Robert Marskar
*/

#include "profile_cylinder_if.H"
#include "cylinder_if.H"

#include <TorusIF.H>
#include <SmoothUnion.H>
#include <SmoothIntersection.H>


profile_cylinder_if::profile_cylinder_if(const RealVect  a_center,
					 const Real      a_length,
					 const Real      a_radius,
					 const int       a_num_left,
					 const int       a_num_right,
					 const Real      a_rad,
					 const Real      a_offset,
					 const Real      a_shift,
					 const Real      a_dist,
					 const Real      a_curv,
					 const bool      a_inside){
#if CH_SPACEDIM==2
  MayDay::Abort("profile_cylinder_if::profile_cylinder_if - this is a 3D geometry only");
#endif

  const RealVect zhat = RealVect(BASISV(2));
  const RealVect c1   = a_center + 0.5*a_length*zhat;
  const RealVect c2   = a_center - 0.5*a_length*zhat;

  Vector<BaseIF*> parts;
  
  BaseIF* cyl = (BaseIF*) (new cylinder_if(c1, c2, a_radius, a_inside));
  parts.push_back(cyl);

  const Real R = a_radius;// - a_rad + a_offset;
  const Real r = a_rad;

  // "Left" profiles
  for (int ileft = 0; ileft < a_num_left; ileft++){
    const RealVect p  = a_center - (ileft+0.5)*(a_dist+2*a_rad)*zhat+ a_shift*zhat;

    BaseIF* torus = (BaseIF*) (new TorusIF(R, r, p, !a_inside));
    parts.push_back(torus);
  }

  // "Right" profiles
  for (int iright = 0; iright < a_num_right; iright++){
    const RealVect p  = a_center + (iright+0.5)*(a_dist+2*a_rad)*zhat+ a_shift*zhat;

    BaseIF* torus = (BaseIF*) (new TorusIF(R, r, p, !a_inside));
    parts.push_back(torus);
  }

  m_baseif = RefCountedPtr<BaseIF> (new SmoothUnion(parts, a_curv));

}

profile_cylinder_if::profile_cylinder_if(const profile_cylinder_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}
  
profile_cylinder_if::~profile_cylinder_if(){

}

Real profile_cylinder_if::value(const RealVect& a_pos) const{
  return m_baseif->value(a_pos);
}

BaseIF* profile_cylinder_if::newImplicitFunction() const{
    return static_cast<BaseIF*> (new profile_cylinder_if(*this));
}
