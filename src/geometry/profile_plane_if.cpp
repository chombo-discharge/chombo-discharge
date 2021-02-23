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
				   const BaseIF*   a_impFunc,
				   const int       a_num_left,
				   const int       a_num_right,
				   const Real      a_ccDist,
				   const Real      a_xShift,
				   const Real      a_yShift,
				   const Real      a_curv,
				   const bool      a_inside){

  const RealVect xhat  = BASISREALV(0);
  const RealVect yhat  = BASISREALV(1);

  const RealVect point = a_point - 1.234E-8*yhat; // Hack to make sure box does not align with grid lines

  // Base planes
#if 0
  BaseIF* base_plane = static_cast<BaseIF*> (new PlaneIF(yhat,  point, a_inside));
  BaseIF* left_plane = static_cast<BaseIF*> (new PlaneIF(-xhat, point - 0.5*a_width*xhat - xhat*1.E-10, a_inside));
  BaseIF* righ_plane = static_cast<BaseIF*> (new PlaneIF(+xhat, point + 0.5*a_width*xhat + xhat*1.E-10, a_inside));
  
  Vector<BaseIF*> parts;
  parts.push_back(base_plane);
  parts.push_back(left_plane);
  parts.push_back(righ_plane);
#else
  Vector<BaseIF*> parts;

  const RealVect lo = a_point - 0.5*a_width*xhat - 1.E3*yhat;
  const RealVect hi = a_point + 0.5*a_width*xhat;
  BaseIF* box = (BaseIF*) (new box_if(lo, hi, false));

  parts.push_back(box);
#endif

  

  // Left profile holes
  for (int ileft = 0; ileft < a_num_left; ileft++){
    TransformIF* transIF = new TransformIF(*a_impFunc);

    const RealVect shift = -(ileft + 0.5)*a_ccDist*xhat + a_xShift*xhat + a_yShift*yhat;
    transIF->translate(shift);

    parts.push_back(transIF);
  }

  // Right profile holes
  for (int ileft = 0; ileft < a_num_left; ileft++){
    TransformIF* transIF = new TransformIF(*a_impFunc);

    const RealVect shift = (ileft + 0.5)*a_ccDist*xhat + a_xShift*xhat + a_yShift*yhat;
    transIF->translate(shift);

    parts.push_back(transIF);
  }

  m_baseif = RefCountedPtr<BaseIF> (new SmoothUnion(parts, a_curv));
}

profile_plane_if::profile_plane_if(const profile_plane_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}

profile_plane_if::~profile_plane_if(){

}


Real profile_plane_if::value(const RealVect& a_pos) const{
  return m_baseif->value(a_pos);
}

BaseIF* profile_plane_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new profile_plane_if(*this));
}
