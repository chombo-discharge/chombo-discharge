#include "Particle.H"
#include "NamespaceHeader.H"

/// default constructor
Particle::Particle() :
  BinItem()
{}

Particle::~Particle()
{}

Particle::Particle(const Real      a_mass,
                   const RealVect& a_position,
                   const RealVect& a_velocity)
  :
  BinItem(a_position),
  m_mass(a_mass),
  m_velocity(a_velocity)
{
  setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
}

void Particle::define(const Real      a_mass,
                      const RealVect& a_position,
                      const RealVect& a_velocity)
{
  setMass(a_mass);
  setPosition(a_position);
  setVelocity(a_velocity);
  setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
}

void Particle::setMass(const Real a_mass)
{
  m_mass = a_mass;
}

const Real& Particle::mass() const
{
  return m_mass;
}

Real& Particle::mass()
{
  return m_mass;
}

// velocity functions
void Particle::setVelocity(const RealVect& a_velocity)
{
  m_velocity = a_velocity;
}

void Particle::setVelocity(const Real& a_velocity,
                            const int   a_dir)
{
  m_velocity[a_dir] = a_velocity;
}

RealVect& Particle::velocity()
{
  return m_velocity;
}

const RealVect& Particle::velocity() const
{
  return m_velocity;
}

Real Particle::velocity(const int a_dir) const
{
  return m_velocity[a_dir];
}

// acceleration functions
void Particle::setAcceleration(const RealVect& a_acceleration)
{
  m_acceleration= a_acceleration;
}

RealVect& Particle::acceleration()
{
  return m_acceleration;
}

const RealVect& Particle::acceleration() const
{
  return m_acceleration;
}

Real Particle::acceleration(const int a_dir) const
{
  return m_acceleration[a_dir];
}

bool Particle::operator == (const Particle& a_p) const
{
  return ( m_mass      == a_p.m_mass     &&
           m_position  == a_p.m_position &&
           m_velocity  == a_p.m_velocity);
}

bool Particle::operator == (const Particle* a_p) const
{
  return (*this == *a_p);
}

bool Particle::operator != (const Particle& a_p) const
{
  return !(*this == a_p);
}

int Particle::size() const
{
  return ( BinItem::size() + sizeof(m_mass) + sizeof(m_velocity) + sizeof(m_acceleration));
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void Particle::linearOut(void* buf) const
{
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

  D_TERM6( *buffer++ = m_acceleration[0];,
	   *buffer++ = m_acceleration[1];,
	   *buffer++ = m_acceleration[2];,
	   *buffer++ = m_acceleration[3];,
	   *buffer++ = m_acceleration[4];,
	   *buffer++ = m_acceleration[5];);

    *buffer = m_mass;
}


// Read a linear (binary) representation of the internal data.
// Assumes the buffer contains the correct data.
void Particle::linearIn(void* buf)
{
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

  D_TERM6( m_acceleration[0] = *buffer++;,
	   m_acceleration[1] = *buffer++;,
	   m_acceleration[2] = *buffer++;,
	   m_acceleration[3] = *buffer++;,
	   m_acceleration[4] = *buffer++;,
	   m_acceleration[5] = *buffer++;);

  m_mass = *buffer;
}

std::ostream & operator<<(std::ostream& ostr, const Particle& p)
{
  ostr << " Particle : " << std::endl;
  ostr << " mass " << p.mass() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.velocity(i); }
  ostr << " ) ";
  ostr << std::endl << " acceleration ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.acceleration(i); }
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
