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
	       const bool& a_inside){
  m_loCorner  = a_loCorner;
  m_hiCorner  = a_hiCorner;
  m_inside    = a_inside;
}

box_if::box_if(const box_if& a_inputIF){
  m_loCorner  = a_inputIF.m_loCorner;
  m_hiCorner  = a_inputIF.m_hiCorner;
  m_inside    = a_inputIF.m_inside;
}


box_if::~box_if(){
  
}

Real box_if::value(const RealVect& a_pos) const{
#if 0 // Original code. Not a signed distance function. 
  RealVect a_point = a_pos;

  a_point -= m_loCorner;
  RealVect newHiCorner = m_hiCorner - m_loCorner;
  Real retval = 1.;
  if(a_point >= RealVect::Zero && a_point <= newHiCorner){
    retval = -1.;
  }

  if(!m_inside) {
    retval = -retval;
  }

  return retval;
#else

  // TLDR: Min(0.0, Max(dx, dy)) is the shortest distance from a point inside the rectangle to one of the edges. It becomes zero if dx > 0 or dy > 0. The
  //       second term max(RealVect::Zero, delta).vectorLength() is funky. In principle we should only take the distance to be sqrt(dx^2 + dy^2 + dz^2) if
  //       we are closest to a "corner", but if we are closest to one of the edges we have dx*dy < 0.0. We can just set the other component to zero, which is what
  //       max(RealVect::Zero, v) does. Then take the distance anyways. 

  const RealVect delta = RealVect(D_DECL(Max(m_loCorner[0] - a_pos[0], a_pos[0] - m_hiCorner[0]),
					 Max(m_loCorner[1] - a_pos[1], a_pos[1] - m_hiCorner[1]),
					 Max(m_loCorner[2] - a_pos[2], a_pos[2] - m_hiCorner[2])));

  const Real retval =  Min(0.0, delta[delta.maxDir(false)]) + max(RealVect::Zero, delta).vectorLength();
  
  return m_inside ? retval : -retval;
#endif


}

BaseIF* box_if::newImplicitFunction() const{
  box_if* boxPtr = new box_if(m_loCorner,m_hiCorner,m_inside);
  return static_cast<BaseIF*> (boxPtr);
}
