/*!
  @file   godunov_particle.cpp
  @brief  Implementation of godunov_particle.H
  @author Robert Marskar
  @date   April 2020
*/

#include "godunov_particle.H"

namespace ChomboDischarge {
  using namespace physics::ito_plasma;

  godunov_particle::godunov_particle() : BinItem(){
  }


  godunov_particle::godunov_particle(const RealVect a_position, const Real a_mass) {
    m_mass     = a_mass;
    m_position = a_position;
  }

  godunov_particle::~godunov_particle(){

  }

  void godunov_particle::define(const RealVect a_position, const Real a_mass){
    setMass(a_mass);
    setPosition(a_position);
  }

  void godunov_particle::setMass(const Real a_mass){
    m_mass = a_mass;
  }

  Real& godunov_particle::mass(){
    return m_mass;
  }

  const Real& godunov_particle::mass() const{
    return m_mass;
  }

  int godunov_particle::size() const{
    return ( BinItem::size() + sizeof(m_mass));
  }

  void godunov_particle::linearOut(void* buf) const{
    Real* buffer = (Real*)buf;
    D_TERM6( *buffer++ = m_position[0];,
	     *buffer++ = m_position[1];,
	     *buffer++ = m_position[2];,
	     *buffer++ = m_position[3];,
	     *buffer++ = m_position[4];,
	     *buffer++ = m_position[5];);

    *buffer = m_mass;
  }

  void godunov_particle::linearIn(void* buf){
    Real* buffer = (Real*)buf;
    D_TERM6( m_position[0] = *buffer++;,
	     m_position[1] = *buffer++;,
	     m_position[2] = *buffer++;,
	     m_position[3] = *buffer++;,
	     m_position[4] = *buffer++;,
	     m_position[5] = *buffer++;);

    m_mass = *buffer;
  }
}
