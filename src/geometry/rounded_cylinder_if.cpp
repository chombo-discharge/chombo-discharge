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

#include "cylinder_if.H"
#include "rounded_cylinder_if.H"

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
  this->makeBaseIF2D();
#elif CH_SPACEDIM==3
  this->makeBaseIF3D();
#endif
}

#if CH_SPACEDIM==2
void rounded_cylinder_if::makeBaseIF2D(){
  MayDay::Abort("rounded_cylinder_if::makeBaseIF2D - not implemented (yet)");
}
#endif

#if CH_SPACEDIM==3
void rounded_cylinder_if::makeBaseIF3D(){

  // TLDR: Construct m_baseif from a main cylinderk.  on each we put a torus and then a smaller cylinder between everything. Default orientation
  //       is along +x, then we rotate. This is constructed such that the baseif function gives a value 

  const RealVect x0 = RealVect::Zero;
  const RealVect x1 = x0 + m_length*BASISREALV(2);
  
  const RealVect y0 = x0 + m_curv*BASISREALV(2);
  const RealVect y1 = x1 - m_curv*BASISREALV(2);

  const Real majorRadius = m_radius - m_curv;
  const Real minorRadius = m_curv;

  BaseIF* mainCylinder   = (BaseIF*) (new cylinder_if(y0, y1, m_radius,        false));
  BaseIF* insideCylinder = (BaseIF*) (new cylinder_if(x0, x1, m_radius-m_curv, false));

  BaseIF* torusBottom = (BaseIF*) (new TorusIF(majorRadius, minorRadius, y0, false));
  BaseIF* torusTop    = (BaseIF*) (new TorusIF(majorRadius, minorRadius, y1, false));

  // Make the intersection of these
  Vector<BaseIF*> parts;
  parts.push_back(mainCylinder);
  parts.push_back(insideCylinder);
  parts.push_back(torusBottom);
  parts.push_back(torusTop);
  BaseIF* isect = (BaseIF*) (new IntersectionIF(parts));

  // Do a transform
  TransformIF* transif = new TransformIF(*isect);

  if(m_center2[2] > m_center1[2]){
    transif->rotate(BASISREALV(2), m_center2-m_center1);
    transif->translate(m_center1);
  }
  else{
    transif->rotate(BASISREALV(2), m_center1-m_center2);
    transif->translate(m_center2);
  }


  // Ok, we're done. 
  m_baseif = RefCountedPtr<BaseIF> (transif);
}
#endif

