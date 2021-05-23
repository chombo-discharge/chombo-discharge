/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointParticle.cpp
  @brief  Implementation of CD_PointParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_PointParticle.H>
#include <CD_NamespaceHeader.H>
  
PointParticle::PointParticle() : BinItem(){
}


PointParticle::PointParticle(const RealVect a_position, const Real a_mass) {
  m_mass     = a_mass;
  m_position = a_position;
}

PointParticle::~PointParticle(){

}

void PointParticle::define(const RealVect a_position, const Real a_mass){
  setMass(a_mass);
  setPosition(a_position);
}

void PointParticle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& PointParticle::mass(){
  return m_mass;
}

const Real& PointParticle::mass() const{
  return m_mass;
}

int PointParticle::size() const{
  return ( BinItem::size() + sizeof(m_mass));
}

void PointParticle::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);

  *buffer = m_mass;
}

void PointParticle::linearIn(void* buf){
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);

  m_mass = *buffer;
}

#include <CD_NamespaceFooter.H>
