/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhysicsLEA.H
  @brief  Interface file for LEA-based ito plasma physics
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaPhysicsLEA_H
#define CD_ItoPlasmaPhysicsLEA_H

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ItoPlasmaPhysicsLEA : public ItoPlasmaPhysics
    {
    public:
      /*!
	@brief Ctor. Does nothing.
      */
      ItoPlasmaPhysicsLEA() {}

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ItoPlasmaPhysicsLEA() {}

      /*!
	@brief Compute Townsend ionization length
      */
      virtual Real
      computeAlpha(const RealVect a_E) const = 0;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      virtual Real
      initialSigma(const Real a_time, const RealVect a_pos) const
      {
        return 0.0;
      }

      /*!
	@brief This must be overridden.
      */
      virtual void
      updateReactionRatesLEA(const RealVect     a_E,
                             const Vector<Real> a_mean_energies,
                             const Real         a_dx,
                             const Real         a_kappa) const override = 0;

      // Disabled in LEA
      Vector<Real>
      computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const override final
      {
        MayDay::Abort("ItoPlasmaPhysicsLEA::computeItoMobilitiesLFA - should not end up here with an LEA model");
      }

      // Disabled in LEA
      Vector<Real>
      computeItoDiffusionLFA(const Real         a_time,
                             const RealVect     a_pos,
                             const RealVect     a_E,
                             const Vector<Real> a_cdr_densities) const override final
      {
        MayDay::Abort("ItoPlasmaPhysicsLEA::compute_diffusion_lfa - should not end up here with an LEA model");
      }

      // Disables in LEA
      void
      updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const override final
      {
        MayDay::Abort("ItoPlasmaPhysicsLEA::updateReactionRatesLFA - should not end up here with an LEA model");
      }
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
