/*!
  @file box_if.cpp
  @brief Contains implementation of the box_if class
  @author Robert Marskar
  @date Nov. 2017
*/

#include "box_if.H"

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
}

BaseIF* box_if::newImplicitFunction() const{
  box_if* boxPtr = new box_if(m_loCorner,m_hiCorner,m_inside);
  return static_cast<BaseIF*> (boxPtr);
}
