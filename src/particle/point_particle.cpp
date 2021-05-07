/*!
  @file   point_particle.cpp
  @brief  Implementation of point_particle.H
  @author Robert Marskar
  @date   April 2020
*/

#include "point_particle.H"

#include "CD_NamespaceHeader.H"
  
point_particle::point_particle() : BinItem(){
}


point_particle::point_particle(const RealVect a_position, const Real a_mass) {
  m_mass     = a_mass;
  m_position = a_position;
}

point_particle::~point_particle(){

}

void point_particle::define(const RealVect a_position, const Real a_mass){
  setMass(a_mass);
  setPosition(a_position);
}

void point_particle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& point_particle::mass(){
  return m_mass;
}

const Real& point_particle::mass() const{
  return m_mass;
}

int point_particle::size() const{
  return ( BinItem::size() + sizeof(m_mass));
}

void point_particle::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);

  *buffer = m_mass;
}

void point_particle::linearIn(void* buf){
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);

  m_mass = *buffer;
}
#include "CD_NamespaceFooter.H"
