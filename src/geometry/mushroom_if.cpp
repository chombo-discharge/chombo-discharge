/*!
  @file    mushroom_if.H
  @brief   Implementation of mushroom_if.H
  @author  Robert Marskar
  @date    Nov. 2017
*/

#include "mushroom_if.H"
#include "rounded_box_if.H"
#include "box_if.H"
#include "cylinder_if.H"
#include "disk3D_if.H"

#include <SmoothIntersection.H>
#include <LatheIF.H>

mushroom_if::mushroom_if(const RealVect a_center,
			 const Real     a_R,
			 const Real     a_r,
			 const Real     a_L,
			 const Real     a_d,
			 const Real     a_curv,
			 const bool     a_inside){

#if CH_SPACEDIM==2
  const RealVect p1 = a_center + a_R*RealVect(BASISV(0));
  const RealVect p2 = a_center - a_R*RealVect(BASISV(0)) - a_d*RealVect(BASISV(1));
  const RealVect p3 = a_center + a_r*RealVect(BASISV(0)) - a_d*RealVect(BASISV(1));
  const RealVect p4 = a_center - a_r*RealVect(BASISV(0)) - (a_d + a_L)*RealVect(BASISV(1));

  BaseIF* box1 = (BaseIF*) (new rounded_box_if(p2, p1, a_curv, !a_inside));
  BaseIF* box2 = (BaseIF*) (new box_if(p4, p3, !a_inside));

  Vector<BaseIF*> parts;
  parts.push_back(box1);
  parts.push_back(box2);

  m_baseif = RefCountedPtr<BaseIF> (new SmoothIntersection(parts, a_curv));

  delete box1;
  delete box2;
  
#elif CH_SPACEDIM==3
  const RealVect zvec = RealVect(0,0,1);

  BaseIF* head = (BaseIF*) (new disk3D_if(a_center, a_R, a_d, a_curv, !a_inside));
  BaseIF* stem = (BaseIF*) (new cylinder_if(a_center - 0.5*a_d*zvec, a_center - (a_d+a_L)*zvec, a_r, !a_inside));

  Vector<BaseIF*> parts;
#if 1 // Original code
  parts.push_back(stem);
  parts.push_back(head);

  m_baseif = RefCountedPtr<BaseIF> (new SmoothIntersection(parts, a_curv));
#else // Debug
  m_baseif = RefCountedPtr<BaseIF> (head);
#endif
#endif
}

mushroom_if::mushroom_if(const mushroom_if& a_inputIF){
  m_baseif = a_inputIF.m_baseif;
}

mushroom_if::~mushroom_if(){

}

Real mushroom_if::value(const RealVect& a_pos) const{
  return m_baseif->value(a_pos);
  MayDay::Abort("stop");
}

BaseIF* mushroom_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new mushroom_if(*this));
}
  
