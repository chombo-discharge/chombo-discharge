/*!
  @file profile_plane_if.H
  @author Robert Marskar
  @date April 2019
*/

#include <IntersectionIF.H>
#include <SphereIF.H>
#include <PlaneIF.H>
#include <TransformIF.H>
#include <UnionIF.H>
#include <SmoothUnion.H>

#include "profile_plane_if.H"
#include "new_sphere_if.H"
#include "box_if.H"

profile_plane_if::profile_plane_if(const RealVect  a_point,
				   const Real      a_width,
				   const int       a_num_left,
				   const int       a_num_right,
				   const Real      a_rad,
				   const Real      a_offset,
				   const Real      a_shift,
				   const Real      a_dist,
				   const Real      a_curv,
				   const bool      a_inside){
  const RealVect xhat = RealVect(BASISV(0));
  const RealVect yhat = RealVect(BASISV(1));

  const RealVect point = a_point - 1.234E-8*yhat; // Hack to make sure box does not align with grid lines
  
  //
  BaseIF* base_plane = static_cast<BaseIF*> (new PlaneIF(yhat,  point, a_inside));
  BaseIF* left_plane = static_cast<BaseIF*> (new PlaneIF(-xhat, point - 0.5*a_width*xhat - xhat*1.E-10, a_inside));
  BaseIF* righ_plane = static_cast<BaseIF*> (new PlaneIF(+xhat, point + 0.5*a_width*xhat + xhat*1.E-10, a_inside));
  
  Vector<BaseIF*> parts;
  parts.push_back(base_plane);
  parts.push_back(left_plane);
  parts.push_back(righ_plane);

  // Left holes
  for (int ileft = 0; ileft < a_num_left; ileft++){
    const RealVect p  = point - (ileft+0.5)*(a_dist+2*a_rad)*xhat + a_offset*yhat + a_shift*xhat;

    BaseIF* sph = (BaseIF*) (new new_sphere_if(p, a_rad, a_inside));
    parts.push_back(sph);
  }

  // Right holes
  for (int iright = 0; iright < a_num_right; iright++){
    const RealVect p  = point + (iright+0.5)*(a_dist+2*a_rad)*xhat + a_offset*yhat + a_shift*xhat;

    BaseIF* sph = (BaseIF*) (new new_sphere_if(p, a_rad, a_inside));
    parts.push_back(sph);
  }

  m_baseif = RefCountedPtr<BaseIF> (new SmoothUnion(parts, a_curv));
}

profile_plane_if::profile_plane_if(const profile_plane_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}

profile_plane_if::~profile_plane_if(){

}


Real profile_plane_if::value(const RealVect& a_pos) const{
  m_baseif->value(a_pos);
}

BaseIF* profile_plane_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new profile_plane_if(*this));
}
