/*!
  @file   ito_particle.H
  @brief  Declaration of a particle class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
  @todo   When we have re-written MeshInterp, we should get rid of the setAcceleration stuff
*/

#ifndef _ITO_PARTICLE_
#define _ITO_PARTICLE_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"

/*!
  @brief A particle class for use with Ito diffusion, i.e. drift-diffusion brownian walkers
*/
class ito_particle : public BinItem {
public:

  /*!
    @brief Default constructor
  */
  ito_particle();

  /*!
    @brief Constructor. Subsequently calls define
  */
  ito_particle(const Real      a_mass,
	       const RealVect& a_position,
	       const RealVect& a_velocity = RealVect::Zero,
	       const Real      a_diffusion = 0.0);

  /*!
    @brief Default dtor
  */
  virtual ~ito_particle();

  /*!
    @brief For non-default ctor
  */
  void define(const Real      a_mass,
              const RealVect& a_position,
              const RealVect& a_velocity = RealVect::Zero,
	      const Real      a_diffusion = 0.0);

  /*!
    @brief Set particle mass
  */
  void setMass(const Real a_mass);

  /*!
    @brief Get mass
  */
  Real& mass();

  /*! 
    @brief Const retrieve mass
  */
  const Real& mass() const;

  /*!
    @brief Set particle diffusion
  */
  void setDiffusion(const Real a_diffusion);

  /*!
    @brief Get diffusion
  */
  Real& diffusion();

  /*! 
    @brief Const retrieve diffusion
  */
  const Real& diffusion() const;

  /*!
    @brief Set the velocity
  */
  void setVelocity(const RealVect& a_velocity);

  /*!
    @brief Set the velocity
  */
  void setVelocity(const Real& a_velocity, const int a_dir);

  /*!
    @brief This is a fucking hack that sets the velocity 
  */
  void setAcceleration(const RealVect& a_velocity){
    setVelocity(a_velocity);
  }

  /*!
    @brief Get reference velocity field
  */
  RealVect& velocity();

  /*!
    @brief Get reference velocity field. Const version. 
  */
  const RealVect& velocity() const;

  /*!
    @brief Get reference velocity field. Const version for one component only. 
  */
  Real velocity(const int a_dir) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const ito_particle& a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator==(const ito_particle* a_p) const;

  /*!
    @brief Comparison operator
  */
  bool operator!=(const ito_particle& a_p) const;

  /*!
    @brief Returns the size, in number of bytes, of a flat representation of the data in this object.
  */
  virtual int size() const;

  /*!
    @brief Write a linear binary representation of the internal data. Assumes that
    sufficient memory for the buffer has already been allocated by
    the caller.
  */
  virtual void linearOut(void* buf) const;

  /*!
    @brief Read a linear binary representation of the internal data.
    Assumes that the buffer has the correct data.
  */
  virtual void linearIn(void* buf);

  /*!
    @brief For outputting particle
  */
  friend std::ostream & operator<<(std::ostream& ostr,const ito_particle& p);

protected:

  Real     m_mass;
  RealVect m_velocity;
  Real     m_diffusion;
};

std::ostream & operator<<(std::ostream& ostr, const ito_particle& p);

#endif

