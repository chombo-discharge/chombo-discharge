/*!
  @file   morrow_bourdon.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _MORROW_BOURDON_
#define _MORROW_BOURDON_

#include "cdr_plasma_physics.H"
#include <CD_PerlinSdf.H>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {
  
    /*!
      @brief Class that implements the Morrow-Lowke discharge model for air.
      @details This class uses the Morrow-Lowke kinetic scheme for air and the three-group Bourdon model for the RTE equations. 
    */
    class morrow_bourdon : public cdr_plasma_physics {
    public:

      morrow_bourdon();
      ~morrow_bourdon();

      void advance_reaction_network(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_Photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_Photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const;


      Vector<Real> compute_cdr_diffusion_coefficients(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const;


      Vector<RealVect> compute_cdr_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const;


      Vector<Real> compute_cdr_electrode_fluxes(const Real         a_time,
						const RealVect     a_pos,
						const RealVect     a_normal,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities,
						const Vector<Real> a_cdr_velocities,
						const Vector<Real> a_cdr_gradients,
						const Vector<Real> a_rte_fluxes,
						const Vector<Real> a_extrap_cdr_fluxes) const;

      Vector<Real> compute_cdr_dielectric_fluxes(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_normal,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities,
						 const Vector<Real> a_cdr_velocities,
						 const Vector<Real> a_cdr_gradients,
						 const Vector<Real> a_rte_fluxes,
						 const Vector<Real> a_extrap_cdr_fluxes) const;

      Vector<Real> compute_cdr_domain_fluxes(const Real           a_time,
					     const RealVect       a_pos,
					     const int            a_dir,
					     const Side::LoHiSide a_side,
					     const RealVect       a_E,
					     const Vector<Real>   a_cdr_densities,
					     const Vector<Real>   a_cdr_velocities,
					     const Vector<Real>   a_cdr_gradients,
					     const Vector<Real>   a_rte_fluxes,
					     const Vector<Real>   a_extrap_cdr_fluxes) const;

      Real initial_sigma(const Real a_time, const RealVect a_pos) const;



      class electron : public CdrSpecies {
      public:

	electron();
	~electron();

	Real initialData(const RealVect a_pos, const Real a_time) const;
    
	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class positive_species : public CdrSpecies {
      public:
	positive_species();
	~positive_species();

	Real initialData(const RealVect a_pos, const Real a_time) const;

	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class negative_species : public CdrSpecies {
      public:
	negative_species();
	~negative_species();

	Real initialData(const RealVect a_pos, const Real a_time) const;

	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class Photon_one : public RtSpecies {
      public:
	Photon_one();
	~Photon_one();
	Real getKappa(const RealVect a_pos) const;
	Real get_lambda() const { return m_lambda;}
	Real get_A()      const { return m_A;}
	Real get_pO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };

      class Photon_two : public RtSpecies {
      public:
	Photon_two();
	~Photon_two();
	Real getKappa(const RealVect a_pos) const;
	Real get_lambda() const { return m_lambda;}
	Real get_A()      const { return m_A;}
	Real get_pO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };

      class Photon_three : public RtSpecies {
      public:
	Photon_three();
	~Photon_three();
	Real getKappa(const RealVect a_pos) const;
	Real get_lambda() const { return m_lambda;}
	Real get_A()      const { return m_A;}
	Real get_pO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };


      Vector<int> m_wallBc;

      int m_nelec_idx; // Electron index
      int m_nplus_idx; // Positive ion index
      int m_nminu_idx; // Negative ion index


      int m_Photon1_idx;
      int m_Photon2_idx;
      int m_Photon3_idx;
      int m_noise_octaves;
      Real m_background_rate;
      Real m_townsend2_conductor;
      Real m_townsend2_dielectric;
      Real m_electrode_yield;
      Real m_dielectric_yield;
      Real m_dielectric_work;
      Real m_N;
      Real m_fracO2;
      Real m_fracN2;
      Real m_pq;
      Real m_p;
      Real m_exc_eff;
      Real m_photo_eff;
      Real m_temp;


      Vector<Real> compute_cathode_flux(const Vector<Real> a_extrapolated_fluxes,
					const Vector<Real> a_ion_densities,
					const Vector<Real> a_ion_velocities,
					const Vector<Real> a_Photon_fluxes,
					const RealVect     a_E,
					const RealVect     a_pos,
					const RealVect     a_normal,
					const Real         a_time) const;

      Vector<Real> compute_anode_flux(const Vector<Real> a_extrapolated_fluxes,
				      const Vector<Real> a_ion_densities,
				      const Vector<Real> a_ion_velocities,
				      const Vector<Real> a_Photon_fluxes,
				      const RealVect     a_E,
				      const RealVect     a_pos,
				      const RealVect     a_normal,
				      const Real         a_time) const;



      Real compute_alpha(const RealVect a_E) const;
      Real compute_eta(const RealVect a_E) const;
      Real compute_eta2(const RealVect a_E) const;
      Real compute_eta3(const RealVect a_E) const;
      Real compute_beta(const RealVect a_E) const;
      Real compute_De(const RealVect a_E) const;
  
      RealVect compute_ve(const RealVect a_E) const;
      RealVect compute_vp(const RealVect a_E) const;
      RealVect compute_vn(const RealVect a_E) const;

      void instantiate_species();
  
      void parse_gas();
      void parse_photoi();
      void parse_see();
      void parse_bc();

    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
