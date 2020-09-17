/*!
  @file   simple_ito_particle.H
  @brief  Implementation of simple_ito_particle.H
  @author Robert Marskar
  @date   Sept. 2020
*/

#include "simple_ito_particle.H"

simple_ito_particle::simple_ito_particle() : BinItem(){
}


simple_ito_particle::simple_ito_particle(const Real a_mass, const RealVect a_position, const Real a_energy) : BinItem(a_position) {
  m_mass   = a_mass;
  m_energy = a_energy;
}

simple_ito_particle::~simple_ito_particle(){

}


void simple_ito_particle::define(const Real a_mass, const RealVect a_position, const Real a_energy){
  setMass(a_mass);
  setPosition(a_position);
  setEnergy(a_energy);
}

void simple_ito_particle::setMass(const Real a_mass){
  m_mass = a_mass;
}

Real& simple_ito_particle::mass(){
  return m_mass;
}

const Real& simple_ito_particle::mass() const{
  return m_mass;
}

void simple_ito_particle::setEnergy(const Real a_energy){
  m_energy = a_energy;
}

Real& simple_ito_particle::energy(){
  return m_energy;
}

const Real& simple_ito_particle::energy() const{
  return m_energy;
}

bool simple_ito_particle::operator==(const simple_ito_particle& a_p) const{
  return ( m_mass      == a_p.m_mass     &&
	   m_energy    == a_p.m_energy   &&
           m_position  == a_p.m_position);
}

bool simple_ito_particle::operator==(const simple_ito_particle* a_p) const{
  return (*this == *a_p);
}

bool simple_ito_particle::operator!=(const simple_ito_particle& a_p) const{
  return !(*this == a_p);
}

int simple_ito_particle::size() const{
  return ( BinItem::size() + sizeof(m_mass) + sizeof(m_energy));
}

void simple_ito_particle::linearOut(void* buf) const{
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

void simple_ito_particle::linearIn(void* buf){
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

std::ostream & operator<<(std::ostream& ostr, const simple_ito_particle& p){
  ostr << " simple_ito_particle : " << std::endl;
  ostr << " mass " << p.mass() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << " energy " << p.energy() << std::endl;
  return ostr;
}
