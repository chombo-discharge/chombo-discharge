/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleImplem.H
  @brief  Implementation of CD_TracerParticle.H
  @author Robert Marskar
*/

#ifndef CD_TracerParticleImplem_H
#define CD_TracerParticleImplem_H

// Our includes
#include <CD_TracerParticle.H>
#include <CD_NamespaceHeader.H>

inline
TracerParticle::TracerParticle() {

  m_position = RealVect::Zero;  
  m_mass     = 0.0;
  m_velocity = RealVect::Zero;
  m_process  = false;
}

inline  
TracerParticle::TracerParticle(const TracerParticle& a_other) {
  this->m_position = a_other.m_position;
  this->m_mass     = a_other.m_mass;
  this->m_velocity = a_other.m_velocity;
  this->m_process  = a_other.m_process;      
}

inline
TracerParticle::~TracerParticle() {
  
}

inline
void TracerParticle::processParticle() {
  m_process = true;
}

inline
void TracerParticle::unprocessParticle() {
  m_process = false;
}

inline
bool& TracerParticle::process() {
  return m_process;
}

inline
const bool& TracerParticle::process() const {
  return m_process;
}

inline
Real& TracerParticle::mass() {
  return m_mass;
}

inline
const Real& TracerParticle::mass() const {
  return m_mass;
}

inline  
RealVect& TracerParticle::velocity() {
  return m_velocity;
}

inline  
const RealVect& TracerParticle::velocity() const {
  return m_velocity;
}

inline
bool TracerParticle::operator==(const TracerParticle& a_p) const{
  return ( this->m_position == a_p.m_position &&
	   this->m_mass     == a_p.m_mass     &&
	   this->m_velocity == a_p.m_velocity &&	   
	   this->m_process  == a_p.m_process);
}

inline
bool TracerParticle::operator!=(const TracerParticle& a_p) const{
  return !(*this == a_p);
}

inline
int TracerParticle::size() const{
  return BinItem::size() + sizeof(m_mass) + sizeof(m_velocity) + sizeof(m_process);
}

inline
void TracerParticle::linearOut(void* buf) const{
  Real* buffer = (Real*) buf;

  // Linearize position onto buffer. 
  *buffer++ = m_position[0];
  *buffer++ = m_position[1];
#if CH_SPACEDIM==3
  *buffer++ = m_position[2];
#endif

  // Linearize mass onto buffer
  *buffer++ = m_mass;

  // Velocity onto buffer
  *buffer++ = m_velocity[0];
  *buffer++ = m_velocity[1];
#if CH_SPACEDIM==3
  *buffer++ = m_velocity[2];
#endif

  // Linearize the processing flag onto the buffer. 
  bool* buf2 = (bool*) buf;
  *buf2++ = m_process;
}

inline
void TracerParticle::linearIn(void* buf) {
  Real* buffer = (Real*) buf;

  m_position[0] = *buffer++;
  m_position[1] = *buffer++;
#if CH_SPACEDIM==3
  m_position[2] = *buffer++;
#endif

  // Put mass onto buffer
  m_mass = *buffer++;

  // Put velocity onto buffer.
  m_velocity[0] = *buffer++;
  m_velocity[1] = *buffer++;
#if CH_SPACEDIM==3
  m_velocity[2] = *buffer++;
#endif

  // Linearize the processing flag onto the buffer. 
  bool* buf2 = (bool*) buf;
  m_process = *buf2++;
}

inline
std::ostream & operator<<(std::ostream& ostr, const TracerParticle& p) {
  ostr << "TracerParticle : \n";

  // Print position. 
  ostr << "\tPosition  = " << p.position() << "\n";
  ostr << "\tMass      = " << p.mass() << "\n";
  ostr << "\tVelocity  = " << p.velocity() << "\n";
  ostr << "\tProcess   = " << p.process() << "\n";      
  
  return ostr;
}

#include <CD_NamespaceFooter.H>

#endif
