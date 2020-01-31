#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _JOINTPARTICLE_H_
#define _JOINTPARTICLE_H_

#include "Particle.H"
#include "NamespaceHeader.H" 

/// Notes: collection of particles, used to keep track of total mass,
/// center of mass position and velocity and number of collected
/// particles. This is needed for AMR calculations, where JointParticles
/// represent the net contribution of all the particles on finer levels.
/// Class members include: mass, position, velocity, number

class JointParticle : public Particle
{
public:
  ///copy constructor and operator= can be compiler defined.

  /// default constructor
  JointParticle();

  /// same as define()
  JointParticle(const Real&     a_mass,
                const RealVect& a_position,
                const RealVect& a_velocity = RealVect::Zero,
                const size_t&   a_numParticles = 0);

  /// default
  virtual ~JointParticle();

  /// for non-default constructor
  void define(const Real&     a_mass,
              const RealVect& a_position,
              const RealVect& a_velocity = RealVect::Zero,
              const size_t&   a_numParticles= 0);

  ///
  virtual void setNumParticles(const size_t& a_numParticles);

  ///
  virtual size_t numParticles() const;

  /// add a particle to the aggregate
  virtual void addParticle(const Particle* const a_p);

  /// add another joint particle to the aggregate
  virtual void addJointParticle(const JointParticle* const a_p);

  /// returns the size, in number of bytes, of a linear representation of the data
  virtual int size() const;

  /// Write a linear binary representation of the internal data.
  virtual void linearOut(void* buf) const;

  /// Read a linear binary representation of the internal data.
  virtual void linearIn(void* buf);

  ///
  virtual void print( void ) const ;

  ///
  friend std::ostream & operator<<(std::ostream& ostr,const JointParticle& p);

  /// The data
protected:
  // number of joint particles
  size_t  m_numParticles;
};


std::ostream & operator<<(std::ostream& ostr, const JointParticle& p);

#include "NamespaceFooter.H"

#endif
