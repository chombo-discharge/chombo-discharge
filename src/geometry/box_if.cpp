/*!
  @file box_if.cpp
  @brief Contains implementation of the box_if class
  @author Robert Marskar
  @date Nov. 2017
*/

#include "box_if.H"

#include <PlaneIF.H>
#include <SmoothUnion.H>
#include <UnionIF.H>
#include <IntersectionIF.H>

box_if::box_if(const RealVect& a_loCorner, 
	       const RealVect& a_hiCorner,
	       const bool&     a_fluidInside){
  m_loCorner    = a_loCorner;
  m_hiCorner    = a_hiCorner;
  m_fluidInside = a_fluidInside;
}

box_if::box_if(const box_if& a_inputIF){
  m_loCorner    = a_inputIF.m_loCorner;
  m_hiCorner    = a_inputIF.m_hiCorner;
  m_fluidInside = a_inputIF.m_fluidInside;
}

box_if::~box_if(){
  
}

Real box_if::value(const RealVect& a_pos) const{
  // TLDR: Min(0.0, Max(dx, dy)) is the shortest distance from a point inside the rectangle to one of the edges. It becomes zero if dx > 0 or dy > 0. The
  //       second term max(RealVect::Zero, delta).vectorLength() is funky. In principle we should only take the distance to be sqrt(dx^2 + dy^2 + dz^2) if
  //       we are closest to a "corner", but if we are closest to one of the edges we have dx*dy < 0.0. We can just set the other component to zero, which is what
  //       max(RealVect::Zero, v) does. Then take the distance anyways. 

  const RealVect delta = RealVect(D_DECL(Max(m_loCorner[0] - a_pos[0], a_pos[0] - m_hiCorner[0]),
					 Max(m_loCorner[1] - a_pos[1], a_pos[1] - m_hiCorner[1]),
					 Max(m_loCorner[2] - a_pos[2], a_pos[2] - m_hiCorner[2])));

  Real retval =  Min(0.0, delta[delta.maxDir(false)]) + max(RealVect::Zero, delta).vectorLength(); // Negative inside box. NOTE: This is RealVect::max and not std::max
  
  if(!m_fluidInside){ // Flip so lsf is positive inside box and negative outside box. 
    retval = -retval;
  }

  return retval;
}

BaseIF* box_if::newImplicitFunction() const{
  return (BaseIF*) (new box_if(m_loCorner, m_hiCorner, m_fluidInside));
}
