/*!
  @file   morrow_zheleznyak.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _MORROW_ZHELEZNYAK_
#define _MORROW_ZHELEZNYAK_

#include <CD_CdrPlasmaPhysics.H>

#include <random>

#include "CD_NamespaceHeader.H"
namespace Physics::cdr_plasma {
  
  /*!
    @brief Class that implements the fluctuating Morrow-Lowke discharge model for air.
  */
  class morrow_zheleznyak : public CdrPlasmaPhysics {
  public:

    morrow_zheleznyak();
    ~morrow_zheleznyak();

    void advanceReactionNetwork(Vector<Real>&          a_particle_sources,
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


    Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities) const;
  
    Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
					    const RealVect     a_pos,
					    const RealVect     a_E,
					    const Vector<Real> a_cdr_densities) const;
  
    Vector<Real> computeCdrDomainFluxes(const Real           a_time,
					   const RealVect       a_pos,
					   const int            a_dir,
					   const Side::LoHiSide a_side,
					   const RealVect       a_E,
					   const Vector<Real>   a_cdr_densities,
					   const Vector<Real>   a_cdr_velocities,
					   const Vector<Real>   a_cdr_gradients,
					   const Vector<Real>   a_rte_fluxes,
					   const Vector<Real>   a_extrap_cdr_fluxes) const;
  
    Vector<Real> computeCdrElectrodeFluxes(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_normal,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities,
					      const Vector<Real> a_cdr_velocities,
					      const Vector<Real> a_cdr_gradients,
					      const Vector<Real> a_rte_fluxes,
					      const Vector<Real> a_extrap_cdr_fluxes) const;

    Vector<Real> computeCdrDielectricFluxes(const Real         a_time,
					       const RealVect     a_pos,
					       const RealVect     a_normal,
					       const RealVect     a_E,
					       const Vector<Real> a_cdr_densities,
					       const Vector<Real> a_cdr_velocities,
					       const Vector<Real> a_cdr_gradients,
					       const Vector<Real> a_rte_fluxes,
					       const Vector<Real> a_extrap_cdr_fluxes) const;

    Vector<Real> compute_dielectric_fluxes(const Vector<Real> a_extrapolated_fluxes,
					   const Vector<Real> a_ion_densities,
					   const Vector<Real> a_ion_velocities,
					   const Vector<Real> a_Photon_fluxes,
					   const RealVect     a_E,
					   const RealVect     a_pos,
					   const RealVect     a_normal,
					   const Real         a_time) const;

    Real initialSigma(const Real a_time, const RealVect a_pos) const;

  
    Vector<int> m_wallBc; // Wall boundary conditions

    bool m_ssa; // Using SSA or not

    bool m_alpha_corr; // Use alpha correction (Soloviev et. al paper)
  
    Real m_seed;
    Real m_cutoff_poisson;
    Real m_poiss_exp_swap;
    std::mt19937_64* m_rng; 
    std::uniform_real_distribution<Real>* m_udist01;
    std::uniform_real_distribution<Real>* m_udist11;
    std::normal_distribution<Real> m_gauss;
    int m_nelec_idx;      // Electron index
    int m_nplus_idx;      // Positive ion index
    int m_nminu_idx;      // Negative ion index
    int m_Photon1_idx;    // UV Photon index

    bool m_extrap_electrode_ebbc;
    bool m_extrap_dielectric_ebbc;

    Real m_background_rate;               // Background ionization
    Real m_townsend2_electrode;           // 2nd Townsend coefficient on electrodes
    Real m_townsend2_dielectric;          // 2nd Townsend coefficient on dielectrics
    Real m_electrode_quantum_efficiency;  // Quantum yield on electrodes
    Real m_dielectric_quantum_efficiency; // Quantum yield on dielectircs
  
    Real m_N;      // Neutral number density
    Real m_p;      // Pressure
    Real m_T;      // Temperature
    Real m_fracO2; // O2 fraction
    Real m_fracN2; // N2 fraction
    Real m_pO2;    //

    Real m_pq;         // Quenching pressure
    Real m_exc_eff;    // Excitation efficiency
    Real m_photoi_eff; // Photoionization efficiency
    Real m_particle_w; // Particle weights

    // Common function for computing CDR fluxes
    Vector<Real> compute_cdr_fluxes(const Real         a_time,
				    const RealVect     a_pos,
				    const RealVect     a_normal,
				    const RealVect     a_E,
				    const Vector<Real> a_cdr_densities,
				    const Vector<Real> a_cdr_velocities,
				    const Vector<Real> a_cdr_gradients,
				    const Vector<Real> a_rte_fluxes,
				    const Vector<Real> a_extrap_cdr_fluxes,
				    const Real         a_townsend2,
				    const Real         a_quantum_efficiency) const;

    // Tau leaping method
    void network_tau(Vector<Real>&          a_particle_sources,
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

    // stochastic simulation algorithm
    void network_ssa(Vector<Real>&          a_particle_sources,
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

    // stochastic simulation algorithm
    void network_rre(Vector<Real>&          a_particle_sources,
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

  
    enum source_comp{
      ssa = 0,
      tau = 1,
      rre = 2
    };

    source_comp m_scomp;

    void parse_gas_params();
    void parse_see();
    void parseDomainBc();
    void parse_ebbc();
    void parse_reaction_settings();
    void parse_initial_particles();
    void parse_alpha_corr();
    void add_uniform_particles(List<Particle>& a_particles);
    void add_gaussian_particles(List<Particle>& a_particles);

    RealVect randomGaussian(const Real a_rad);

    RealVect randomDirection();
#if CH_SPACEDIM==2
    RealVect randomDirection2D();
#elif CH_SPACEDIM==3
    RealVect randomDirection3D();
#endif
  
    int poisson_reaction(const Real a_propensity, const Real a_dt) const;
  
    Real computeAlpha(const RealVect a_E) const;
    Real compute_eta(const RealVect a_E)   const;
    Real compute_eta2(const RealVect a_E)  const;
    Real compute_eta3(const RealVect a_E)  const;
    Real compute_beta(const RealVect a_E)  const;

    Real compute_De(const RealVect a_E) const;

    RealVect compute_ve(const RealVect a_E) const;
    RealVect compute_vp(const RealVect a_E) const;
    RealVect compute_vn(const RealVect a_E) const;

    Real excitation_rates(const Real a_E) const;
    Real sergey_factor(const Real a_O2frac) const;

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

    class uv_Photon : public RtSpecies {
    public:
      uv_Photon();
      ~uv_Photon();


      Real getKappa(const RealVect a_pos) const;
      Real get_random_kappa() const;

      Real get_pO2() const { return m_pO2; }

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

  };

}
#include "CD_NamespaceFooter.H"
#endif
