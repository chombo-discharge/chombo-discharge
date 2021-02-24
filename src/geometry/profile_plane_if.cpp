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
#include <ComplementIF.H>
#include <SmoothUnion.H>
#include <SmoothIntersection.H>

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
				   const bool      a_fluidInside){

  const RealVect xhat  = BASISREALV(0);
  const RealVect yhat  = BASISREALV(1);

  const RealVect point = a_point - 1.234E-8*yhat; // Hack to make sure box does not align with grid lines

  // Make base box
  Vector<BaseIF*> parts;
  const RealVect lo = point - 0.5*a_width*xhat - 1.E3*yhat;
  const RealVect hi = point + 0.5*a_width*xhat;
  BaseIF* box = (BaseIF*) (new box_if(lo, hi, false)); // Construct base box with fluid outside. 
  //  parts.push_back(box);

  // Left profile holes
  for (int ileft = 0; ileft < a_num_left; ileft++){
    TransformIF* transIF = new TransformIF(*a_impFunc);

    const RealVect shift = -(ileft + 0.5)*a_ccDist*xhat + a_xShift*xhat + a_yShift*yhat;
    transIF->translate(shift);

    parts.push_back(transIF);
  }

  // Right profile holes
  for (int iright = 0; iright < a_num_right; iright++){
    TransformIF* transIF = new TransformIF(*a_impFunc);

    const RealVect shift = (iright + 0.5)*a_ccDist*xhat + a_xShift*xhat + a_yShift*yhat;
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
  Real retval = m_baseif->value(a_pos); 

  if(m_fluidInside){ // m_baseif was constructed such that the fluid was outside. Revert if fluid is inside. 
    retval = -retval;
  }

  return retval;
}

BaseIF* profile_plane_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new profile_plane_if(*this));
}
