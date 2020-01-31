#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"
#include "NamespaceHeader.H"

/// A basic particle class. 
/// Includes mass, position, velocity, and acceleration

class Particle : public BinItem
{
public:

  /// default constructor
  Particle();

  /// same as define()
  Particle(const Real      a_mass,
           const RealVect& a_position,
           const RealVect& a_velocity = RealVect::Zero);

  /// default
  virtual ~Particle();

  /// for non-default constructor
  void define(const Real      a_mass,
              const RealVect& a_position,
              const RealVect& a_velocity = RealVect::Zero);

  /// set particle mass
  void setMass(const Real a_mass);

  // retrieve reference mass field
  Real& mass();

  // retrieve const. reference mass field
  const Real& mass() const;

  /// set velocity field
  void setVelocity(const RealVect& a_velocity);
  void setVelocity(const Real& a_velocity, const int a_dir);

  // retrieve reference velocity field
  RealVect& velocity();

  // retrieve const. reference velocity field
  const RealVect& velocity() const;

  // retrieve vel. comp. in a given direction
  Real velocity(const int a_dir) const;

  /// set acceleration field
  void setAcceleration(const RealVect& a_velocity);

  // retrieve reference acceleration field
  RealVect& acceleration();

  // retrieve const. reference acceleration field
  const RealVect& acceleration() const;

  // retrieve acc. comp. in a given direction
  Real acceleration(const int a_dir) const;

  virtual bool operator==(const Particle& a_p) const;
  virtual bool operator==(const Particle* a_p) const;
  virtual bool operator!=(const Particle& a_p) const;

  /// {\bf linearization functions}

  /** Returns the size, in number of bytes, of a flat linear
      representation of the data in this object.
   */
  virtual int size() const;

  // return size, in number of bytes, of a comp of the internal data
  //virtual size_t sizeOfComp(const int a_comp) const;

  ///
  /** Write a linear binary representation of the internal data. Assumes that
      sufficient memory for the buffer has already been allocated by
      the caller.
  */
  virtual void linearOut(void* buf) const;

  ///
  /** Read a linear binary representation of the internal data.
      Assumes that the buffer has the correct data.
  */
  virtual void linearIn(void* buf);

  friend std::ostream & operator<<(std::ostream& ostr,const Particle& p);

protected:

  Real     m_mass;
  RealVect m_velocity;
  RealVect m_acceleration;
};

std::ostream & operator<<(std::ostream& ostr, const Particle& p);

#include "NamespaceFooter.H"

#endif // include guard

