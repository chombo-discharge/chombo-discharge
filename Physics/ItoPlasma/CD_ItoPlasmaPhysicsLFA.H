/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhysicsLFA.H
  @brief  Interface file for LFA-based ito plasma physics.
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaPhysicsLFA_H
#define CD_ItoPlasmaPhysicsLFA_H

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Middle layer for describing ItoPlasmaPhysics in the LFA
    */
    class ItoPlasmaPhysicsLFA : public ItoPlasmaPhysics {
    public:
      
      /*!
	@brief Ctor. Does nothing.
      */
      ItoPlasmaPhysicsLFA(){

      }

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ItoPlasmaPhysicsLFA(){

      }

      /*!
	@brief Compute Townsend ionization length
      */
      virtual Real computeAlpha(const RealVect a_E) const = 0;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      virtual Real initialSigma(const Real a_time, const RealVect a_pos) const {return 0.0;}

      /*!
	@brief This must be overridden.
      */
      void updateReactionRatesLEA(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const override final {
	MayDay::Abort("ItoPlasmaPhysicsLEA::updateReactionRatesLEA - should not end up here with an LFA model");
      }

      /*!
	@brief Compute mobilities. 
      */
      virtual Vector<Real> computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const override  = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real> computeItoDiffusionLFA(const Real         a_time,
						  const RealVect     a_pos,
						  const RealVect     a_E,
						  const Vector<Real> a_cdr_densities) const override = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const override = 0;


      
    };
  }
}

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaPhysicsImplem.H>

#endif
