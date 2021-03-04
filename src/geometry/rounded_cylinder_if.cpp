/*!
  @file   rounded_cylinder_if.cpp
  @brief  Implementation of rounded_cylinder_if.H
  @date   Feb. 2021
  @author Robert Marskar
*/

#include <PolyGeom.H>
#include <TorusIF.H>
#include <IntersectionIF.H>
#include <TransformIF.H>
#include <PlaneIF.H>
#include <SmoothIntersection.H>
#include <SmoothUnion.H>

#include "cylinder_if.H"
#include "rounded_cylinder_if.H"
#include "rounded_box_if.H"
#include "torus_if.H"

rounded_cylinder_if::rounded_cylinder_if(const RealVect a_center1, const RealVect a_center2, const Real a_radius, const Real a_curv, const bool a_fluidInside){
  m_center1      = a_center1;
  m_center2      = a_center2;
  m_length       = (m_center2-m_center1).vectorLength();
  m_radius       = a_radius;
  m_curv         = a_curv;
  m_fluidInside  = a_fluidInside;

  this->makeBaseIF();
}

rounded_cylinder_if::rounded_cylinder_if(const rounded_cylinder_if& a_inputIF){
  m_fluidInside = a_inputIF.m_fluidInside;
  m_baseif      = a_inputIF.m_baseif;
}

Real rounded_cylinder_if::value(const RealVect& a_point) const {
  Real retval = m_baseif->value(a_point);

  if(m_fluidInside){
    retval = -retval;
  }
  
  return retval;
}

BaseIF* rounded_cylinder_if::newImplicitFunction() const{
  return (BaseIF*) (new rounded_cylinder_if(*this));
}

void rounded_cylinder_if::makeBaseIF(){
#if CH_SPACEDIM==2
  BaseIF* bif = this->makeBaseIF2D();
#elif CH_SPACEDIM==3
  BaseIF* bif = this->makeBaseIF3D();
#endif

  // Rotate and translate into place
  TransformIF* transif = new TransformIF(*bif);

  const int dir     = SpaceDim - 1;
  const RealVect up = BASISREALV(dir);
  if(m_center2[dir] >= m_center1[dir]){
    transif->rotate(up, m_center2-m_center1);
    transif->translate(m_center1);
  }
  else{
    transif->rotate(up, m_center1-m_center2);
    transif->translate(m_center2);
  }

  m_baseif = RefCountedPtr<BaseIF> (transif);
}

#if CH_SPACEDIM==2
BaseIF* rounded_cylinder_if::makeBaseIF2D(){
  const RealVect x0 = RealVect::Zero - m_radius*BASISREALV(0);
  const RealVect x1 = RealVect::Zero + m_radius*BASISREALV(0) + m_length*BASISREALV(1);

  return (BaseIF*) (new rounded_box_if(x0, x1, m_curv, false));
}
#endif

#if CH_SPACEDIM==3
BaseIF* rounded_cylinder_if::makeBaseIF3D(){

  // TLDR: Construct m_baseif from a main cylinderk.  on each we put a torus and then a smaller cylinder between everything. Default orientation
  //       is along +z.

  const RealVect up = BASISREALV(SpaceDim-1);
  
  const RealVect x0 = RealVect::Zero;
  const RealVect x1 = x0 + m_length*up;
  
  const RealVect y0 = x0 + m_curv*up;
  const RealVect y1 = x1 - m_curv*up;

  const Real majorRadius = m_radius - m_curv;
  const Real minorRadius = m_curv;

  BaseIF* mainCylinder   = (BaseIF*) (new cylinder_if(y0, y1, m_radius,        false));
  BaseIF* insideCylinder = (BaseIF*) (new cylinder_if(x0, x1, m_radius-m_curv, false));

  BaseIF* torusBottom = (BaseIF*) (new torus_if(y0, majorRadius, minorRadius, false));
  BaseIF* torusTop    = (BaseIF*) (new torus_if(y1, majorRadius, minorRadius, false));

  // Make the intersection of these
  Vector<BaseIF*> parts;
  parts.push_back(mainCylinder);
  parts.push_back(insideCylinder);
  parts.push_back(torusBottom);
  parts.push_back(torusTop);
  
  return (BaseIF*) (new IntersectionIF(parts));
}
#endif

