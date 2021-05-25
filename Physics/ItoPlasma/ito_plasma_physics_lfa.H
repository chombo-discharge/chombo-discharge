/*!
  @file   ito_plasma_physics_lea.H
  @brief  Interface file for LEA-based ito plasma physics
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_PHYSICS_LFA_
#define _ITO_PLASMA_PHYSICS_LFA_

#include "ito_plasma_physics.H"

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace ito_plasma {

    /*!
      @brief Middle layer for describing ito_plasma_physics in the LFA
    */
    class ito_plasma_physics_lfa : public ito_plasma_physics {
    public:
      
      /*!
	@brief Ctor. Does nothing.
      */
      ito_plasma_physics_lfa(){

      }

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ito_plasma_physics_lfa(){

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
      void update_reaction_rates_lea(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const override final {
	MayDay::Abort("ito_plasma_physics_lea::update_reaction_rates_lea - should not end up here with an LFA model");
      }

      /*!
	@brief Compute mobilities. 
      */
      virtual Vector<Real> compute_ito_mobilities_lfa(const Real a_time, const RealVect a_pos, const RealVect a_E) const override  = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real> compute_ito_diffusion_lfa(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_E,
						     const Vector<Real> a_cdr_densities) const override = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void update_reaction_rates_lfa(const RealVect a_E, const Real a_dx, const Real a_kappa) const override = 0;


      
    };
  }
}
#include "CD_NamespaceFooter.H"

#include "ito_plasma_physicsI.H"

#endif
