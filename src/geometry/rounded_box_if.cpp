/*!
  @file rounded_box_if.cpp
  @brief Implementation of rounded_box_if.H
  @date Nov. 2017
  @author Robert Marskar
*/

#include "rounded_box_if.H"

#include <PlaneIF.H>
#include <SmoothUnion.H>
#include <TransformIF.H>

rounded_box_if::rounded_box_if(const RealVect a_loCorner,
			       const RealVect a_hiCorner,
			       const Real     a_curv,
			       const bool     a_fluidInside){
  m_fluidInside = a_fluidInside;

  const RealVect xyz = a_hiCorner - a_loCorner;
  
  // Make a slab whose center is at zero and whose widths are given by a_xyz
  Vector<BaseIF*> parts;
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int s      = sign(sit());
      const RealVect n = s*BASISREALV(dir);
      const RealVect p = n*0.5*xyz[dir];

      BaseIF* baseif = (BaseIF*) new PlaneIF(n, p, true);
      parts.push_back(baseif);
    }
  }

  // Do rounded corners. 
  BaseIF* bif = (BaseIF*) new SmoothUnion(parts, a_curv);

  // Translate box
  TransformIF* tif = new TransformIF(*bif);
  tif->translate(a_loCorner + 0.5*xyz);

  // Done;
  m_baseif = RefCountedPtr<BaseIF> (tif);
}

rounded_box_if::rounded_box_if(const rounded_box_if& a_inputIF){
  m_fluidInside = a_inputIF.m_fluidInside;
  m_baseif      = a_inputIF.m_baseif;
}

rounded_box_if::~rounded_box_if(){
}

Real rounded_box_if::value(const RealVect& a_point) const{

  // TLDR: m_baseif designed so that f < 0 outside the box. This means that the fluid is outside. Revert if inside. 
  Real retval = m_baseif->value(a_point);

  if(m_fluidInside){
    retval = -retval;
  }

  return retval;
}

BaseIF* rounded_box_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new rounded_box_if(*this));
}

