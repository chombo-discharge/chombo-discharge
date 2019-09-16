/*!
  @file disk3D_if.H
*/

#include "disk3D_if.H"
#include "cylinder_if.H"

#include <TorusIF.H>
#include <SmoothUnion.H>
#include <SmoothIntersection.H>
#include <IntersectionIF.H>

disk3D_if::disk3D_if(const RealVect a_point, const Real a_R, const Real a_d, const Real a_curv, const bool a_inside){
#if CH_SPACEDIM==3
  const RealVect zvec = RealVect(0,0,1);


  BaseIF* cyl1 = (BaseIF*) (new cylinder_if(a_point, a_point - a_d*zvec, a_R - a_curv, a_inside));
  BaseIF* cyl2 = (BaseIF*) (new cylinder_if(a_point - a_curv*zvec, a_point - a_d*zvec + a_curv*zvec, a_R, a_inside));
  BaseIF* tor1 = (BaseIF*) (new TorusIF(a_R-a_curv, a_curv, a_point - a_curv*zvec, a_inside));
  BaseIF* tor2 = (BaseIF*) (new TorusIF(a_R-a_curv, a_curv, a_point - a_d*zvec + a_curv*zvec, a_inside));

  Vector<BaseIF*> parts;
  parts.push_back(cyl1);
  parts.push_back(cyl2);
  parts.push_back(tor1);
  parts.push_back(tor2);
 
  
  m_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts));
#else
  MayDay::Abort("disk3D_if::disk3D_if - this is a 3D class only!");
#endif
}

disk3D_if::~disk3D_if(){

}

disk3D_if::disk3D_if(const disk3D_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}

Real disk3D_if::value(const RealVect& a_pos) const{
  return m_baseif->value(a_pos);
}

BaseIF* disk3D_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new disk3D_if(*this));
}
