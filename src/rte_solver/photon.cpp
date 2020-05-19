/*!
  @file   photon.H
  @brief  Declaration of a photon class for particle methods
  @author Robert Marskar
  @date   May 2019
*/

#include "photon.H"

photon::photon() : BinItem(){

}

photon::~photon(){

}

photon::photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass){
  this->define(a_position, a_velocity, a_kappa, a_mass);
}

void photon::define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass){
  setPosition(a_position);
  setVelocity(a_velocity);
  setKappa(a_kappa);
  setMass(a_mass);
}

// Set get functions
void photon::setKappa(const Real a_kappa)           { m_kappa    = a_kappa;    }
void photon::setMass(const Real a_mass)             { m_mass   = a_mass;   }
void photon::setVelocity(const RealVect& a_velocity){ m_velocity = a_velocity; }

Real&     photon::mass()    { return m_mass;   }
Real&     photon::kappa()   { return m_kappa;    }
RealVect& photon::velocity(){ return m_velocity; }

const Real&     photon::kappa()    const{ return m_kappa;    }
const Real&     photon::mass()     const{ return m_mass;   }
const RealVect& photon::velocity() const{ return m_velocity; }

// Comparison functions
bool photon::operator == (const photon* a_p) const{ return (*this == *a_p); }
bool photon::operator != (const photon& a_p) const{ return !(*this == a_p); }
bool photon::operator == (const photon& a_p) const {
  return ( m_position  == a_p.m_position &&
           m_velocity  == a_p.m_velocity &&
	   m_kappa     == a_p.m_kappa &&
	   m_mass    == a_p.m_mass);
}

int photon::size() const{ return ( BinItem::size() + sizeof(m_velocity) + sizeof(m_kappa) + sizeof(m_mass)); }

void photon::linearOut(void* buf) const{
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);

  D_TERM6( *buffer++ = m_velocity[0];,
	   *buffer++ = m_velocity[1];,
	   *buffer++ = m_velocity[2];,
	   *buffer++ = m_velocity[3];,
	   *buffer++ = m_velocity[4];,
	   *buffer++ = m_velocity[5];);
  
  *buffer = m_kappa;
  buffer++;
  *buffer = m_mass;
}

void photon::linearIn(void* buf){
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);

  D_TERM6( m_velocity[0] = *buffer++;,
	   m_velocity[1] = *buffer++;,
	   m_velocity[2] = *buffer++;,
	   m_velocity[3] = *buffer++;,
	   m_velocity[4] = *buffer++;,
	   m_velocity[5] = *buffer++;);

  m_kappa = *buffer;
  buffer++;
  m_mass = *buffer;
}
