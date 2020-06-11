/*!
  @file   ito_plasma_air2.H
  @brief  Simple class for 
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_AIR2_
#define _ITO_PLASMA_AIR2_

#include "ito_plasma_physics.H"

#include <chrono>
#include <random>

namespace physics {
  namespace ito_plasma{

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_air2 : public ito_plasma_physics {
    public:

      class electron;
      class positive;
      class photonZ;
      
      ito_plasma_air2();
      ~ito_plasma_air2();

      Real compute_alpha(const RealVect a_E) const override;

      Vector<RealVect> compute_ito_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const override;

      Vector<Real> compute_ito_diffusion(const Real         a_time,
					 const RealVect     a_pos,
					 const RealVect     a_E,
					 const Vector<Real> a_cdr_densities) const override;

      void advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
				    Vector<List<photon>* >&       a_photons,
				    Vector<List<photon>* >&       a_newPhotons,
				    const RealVect                a_E,           
				    const RealVect                a_pos,
				    const Real                    a_dx,
				    const Real                    a_kappa, 
				    const Real                    a_dt) const override;

    protected:

      // Stuff for initial conditions
      Real m_blob_radius;
      int m_num_particles;
      int m_seed;
      RealVect m_blob_center;
      mutable std::mt19937_64 m_rng;
      mutable std::uniform_real_distribution<Real> m_udist11;
      mutable std::normal_distribution<Real> m_gauss;

      // Indices for species
      int m_electron_idx;
      int m_positive_idx;
      int m_photonZ_idx;

      // Pressure, neutral number density etc
      Real m_N;      
      Real m_p;      
      Real m_pq;     
      Real m_T;      
      Real m_O2frac; 
      Real m_N2frac; 
      Real m_pO2;    
      Real m_photoi_factor; 

      // Functions for drawing initial particles
      void draw_initial_particles();
      RealVect random_gaussian();
      RealVect random_direction() const;

      Real excitation_rates(const Real a_E) const;
      Real sergey_factor(const Real a_O2frac) const;
    };

    
    class ito_plasma_air2::electron : public ito_species {
    public:
      electron();
      ~electron();
    };

    class ito_plasma_air2::positive : public ito_species {
    public:
      positive();
      ~positive();
    };

    class ito_plasma_air2::photonZ : public rte_species {
    public:
      photonZ();
      ~photonZ();

      Real get_kappa(const RealVect a_pos) const override;
      
    protected:
      Real m_pO2;
      Real m_f1;
      Real m_f2;
      Real m_K1;
      Real m_K2;

      Real m_seed;
      std::mt19937_64* m_rng; 
      std::uniform_real_distribution<Real>* m_udist01;
    };
  }
}
#endif
