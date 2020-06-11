/*!
  @file   ito_plasma_physics.H
  @brief  Main file for describing Ito-based plasma physics
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_PHYSICS_
#define _ITO_PLASMA_PHYSICS_

#include <RealVect.H>
#include <RefCountedPtr.H>
#include <LoHiSide.H>
#include <List.H>

#include "ito_species.H"
#include "rte_species.H"
#include "photon.H"
#include "ito_particle.H"

namespace physics {
  namespace ito_plasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_physics {
    public:
      ito_plasma_physics(){}
      ~ito_plasma_physics(){}

      /*!
	@brief Compute Townsend ionization length
      */
      virtual Real compute_alpha(const RealVect a_E) const = 0;

      /*!
	@brief Get all species
      */
      const Vector<RefCountedPtr<ito_species> >& get_ito_species() const { 
	return m_ito_species; 
      }

      /*!
	@brief Get all photons
      */
      const Vector<RefCountedPtr<rte_species> >& get_rte_species() const {
	return m_rte_species;
      }

      /*!
	@brief Return number of ion equations
      */
      int get_num_ito_species() const{
	return m_ito_species.size();
      }

      /*!
	@brief Return number of RTE equations
      */
      int get_num_rte_species() const {
	return m_rte_species.size();
      }

      /*!
	@brief Set initial surface charge. Override if you want.
      */
      virtual Real initial_sigma(const Real a_time, const RealVect a_pos) const {
	return 0.0;
      }

      /*!
	@brief Compute ito velocities
      */
      virtual Vector<RealVect> compute_ito_velocities(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real> compute_ito_diffusion(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Reaction network advance
      */
      virtual void advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
					    Vector<List<photon>* >&       a_photons,
					    Vector<List<photon>* >&       a_newPhotons,
					    const RealVect                a_E,           
					    const RealVect                a_pos,
					    const Real                    a_dx,
					    const Real                    a_kappa, 
					    const Real                    a_dt) const = 0;

    protected:

      Vector<RefCountedPtr<ito_species> > m_ito_species;
      Vector<RefCountedPtr<rte_species> > m_rte_species;

      int m_num_ito_species;
      int m_num_rte_species;
    };
  }
}
#endif
