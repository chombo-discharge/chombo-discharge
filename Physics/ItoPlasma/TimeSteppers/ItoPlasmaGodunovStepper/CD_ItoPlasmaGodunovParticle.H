/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaGodunovParticle.H
  @brief  Declaration of a particle class for ItoPlasmaGodunovStepper
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaGodunovParticle_H
#define CD_ItoPlasmaGodunovParticle_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Chombo includes
#include <BinItem.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief A particle class for use with Ito diffusion, i.e. drift-diffusion brownian walkers
    */
    class ItoPlasmaGodunovParticle : public BinItem
    {
    public:
      /*!
	@brief Default constructor
      */
      ItoPlasmaGodunovParticle();

      /*!
	@brief Constructor. Subsequently calls define
      */
      ItoPlasmaGodunovParticle(const RealVect a_position, const Real a_mass);

      /*!
	@brief Default dtor
      */
      virtual ~ItoPlasmaGodunovParticle();

      /*!
	@brief For non-default ctor
      */
      void
      define(const RealVect a_position, const Real a_mass);

      /*!
	@brief Set particle mass
      */
      void
      setMass(const Real a_mass);

      /*!
	@brief Get mass
      */
      Real&
      mass();

      /*! 
	@brief Const retrieve mass
      */
      const Real&
      mass() const;

      /*!
	@brief Returns the size, in number of bytes, of a flat representation of the data in this object.
      */
      virtual int
      size() const;

      /*!
	@brief Write a linear binary representation of the internal data. Assumes that
	sufficient memory for the buffer has already been allocated by
	the caller.
      */
      virtual void
      linearOut(void* buf) const;

      /*!
	@brief Read a linear binary representation of the internal data.
	Assumes that the buffer has the correct data.
      */
      virtual void
      linearIn(void* buf);

    protected:
      Real m_mass;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
