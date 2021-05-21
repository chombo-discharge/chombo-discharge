/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_SimpleItoParticle.H
  @brief  Implementation of CD_SimpleItoParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_SimpleItoParticle.H>
#include <CD_NamespaceHeader.H>

SimpleItoParticle::SimpleItoParticle() : BinItem(){
}

SimpleItoParticle::SimpleItoParticle(const Real a_mass, const RealVect a_position, const Real a_energy) : BinItem(a_position) {
  m_mass   = a_mass;
  m_energy = a_energy;
}

SimpleItoParticle::~SimpleItoParticle(){

}

void SimpleItoParticle::define(const Real a_mass, const RealVect a_position, const Real a_energy){
  setMass(a_mass);
  setPosition(a_position);
  setEnergy(a_energy);
}

void SimpleItoParticle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& SimpleItoParticle::mass(){
  return m_mass;
}

const Real& SimpleItoParticle::mass() const{
  return m_mass;
}

void SimpleItoParticle::setEnergy(const Real a_energy){
  m_energy = a_energy;
}

Real& SimpleItoParticle::energy(){
  return m_energy;
}

const Real& SimpleItoParticle::energy() const{
  return m_energy;
}

bool SimpleItoParticle::operator==(const SimpleItoParticle& a_p) const{
  return ( m_mass      == a_p.m_mass     &&
	   m_energy    == a_p.m_energy   &&
	   m_position  == a_p.m_position);
}

bool SimpleItoParticle::operator==(const SimpleItoParticle* a_p) const{
  return (*this == *a_p);
}

bool SimpleItoParticle::operator!=(const SimpleItoParticle& a_p) const{
  return !(*this == a_p);
}

int SimpleItoParticle::size() const{
  return ( BinItem::size() + sizeof(m_mass) + sizeof(m_energy));
}

void SimpleItoParticle::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);

  *buffer++ = m_mass;
  *buffer++ = m_energy;
}

void SimpleItoParticle::linearIn(void* buf){
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);

  m_mass   = *buffer++;
  m_energy = *buffer++;
}

std::ostream & operator<<(std::ostream& ostr, const SimpleItoParticle& p){
  ostr << " SimpleItoParticle : " << std::endl;
  ostr << " mass " << p.mass() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << " energy " << p.energy() << std::endl;
  return ostr;
}

#include <CD_NamespaceFooter.H>
