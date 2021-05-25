/*!
  @file   ito_plasma_physics_lea.H
  @brief  Interface file for LEA-based ito plasma physics
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_PHYSICS_LEA_
#define _ITO_PLASMA_PHYSICS_LEA_

#include <CD_ItoPlasmaPhysics.H>

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_physics_lea : public ItoPlasmaPhysics {
    public:

      /*!
	@brief Ctor. Does nothing.
      */
      ito_plasma_physics_lea(){

      }

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ito_plasma_physics_lea(){

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
      virtual void updateReactionRatesLEA(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const override = 0;

      // Disabled in LEA
      Vector<Real> computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const override final {
	MayDay::Abort("ito_plasma_physics_lea::computeItoMobilitiesLFA - should not end up here with an LEA model");
      }

      // Disabled in LEA
      Vector<Real> computeItoDiffusionLFA(const Real         a_time,
					     const RealVect     a_pos,
					     const RealVect     a_E,
					     const Vector<Real> a_cdr_densities) const override final {
	MayDay::Abort("ito_plasma_physics_lea::compute_diffusion_lfa - should not end up here with an LEA model");
      }

      // Disables in LEA
      void updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const override final {
	MayDay::Abort("ito_plasma_physics_lea::updateReactionRatesLFA - should not end up here with an LEA model");
      }


      
    };
  }
}
#include "CD_NamespaceFooter.H"

#endif
