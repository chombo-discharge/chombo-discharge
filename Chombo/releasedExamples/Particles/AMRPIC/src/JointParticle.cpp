#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cmath>
#include "JointParticle.H"
#include "RealVect.H"
#include "NamespaceHeader.H"

// default constructor
JointParticle::JointParticle()
  :
  Particle(),
  m_numParticles(0)
{}

JointParticle::~JointParticle()
{}

JointParticle::JointParticle(const Real&     a_mass,
                             const RealVect& a_position,
                             const RealVect& a_velocity,
                             const size_t&   a_numParticles)
  :
  Particle(a_mass,a_position,a_velocity)
{
  setNumParticles(a_numParticles);
}

void JointParticle::define(const Real&     a_mass,
                           const RealVect& a_position,
                           const RealVect& a_velocity,
                           const size_t&   a_numParticles)
{
  setMass(a_mass);
  setPosition(a_position);
  setVelocity(a_velocity);
  setNumParticles(a_numParticles);
}

void JointParticle::setNumParticles(const size_t& a_numParticles)
{
  //CH_assert(a_numParticles>=0);
  m_numParticles = a_numParticles;
}

size_t JointParticle::numParticles() const
{
  return m_numParticles;
}

/// add a particle to the aggregate
void JointParticle::addParticle(const Particle* const a_p)
{
  m_position = m_position*m_mass + a_p->position()*a_p->mass();
  m_velocity = m_velocity*m_mass + a_p->velocity()*a_p->mass();
  m_mass    += a_p->mass();

  m_position /= m_mass;
  m_velocity /= m_mass;

  m_numParticles++;
}

/// add another joint particle to the aggregate
void JointParticle::addJointParticle(const JointParticle* const a_p)
{
  m_position = m_position*m_mass + a_p->position()*a_p->mass();
  m_velocity = m_velocity*m_mass + a_p->velocity()*a_p->mass();
  m_mass    += a_p->mass();

  m_position /= m_mass;
  m_velocity /= m_mass;

  m_numParticles += a_p->numParticles();
}


int JointParticle::size() const
{
  return ( Particle::size() + sizeof(size_t) );
}


// diabled cause this class is not meant to be used with this function
void JointParticle::linearOut(void* buf) const
{
  char* buffer = (char*)buf;
  Particle::linearOut( buffer );
  buffer += Particle::size();

  *((size_t*)buffer) = m_numParticles;
}

// diabled cause this class is not meant to be used with this function
void JointParticle::linearIn(void* buf)
{
  char* buffer = (char*)buf;
  Particle::linearIn( buf );
  buffer += Particle::size();

  m_numParticles = *((size_t*)buffer);
}

// Write a human-readable version to stdout.  This is for debugging.
void JointParticle::print(void) const
{
  std::cout << *this ;
}

std::ostream & operator<<(std::ostream& ostr, const JointParticle& p)
{
  ostr << " JointParticle : ";
  ostr << " mass " << p.mass() << std::endl;
  ostr << "position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.velocity(i); }
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
