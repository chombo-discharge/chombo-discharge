/*!
  @file   ito_plasma_air3.H
  @brief  Simple class for 3-species air plasma
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_AIR3_
#define _ITO_PLASMA_AIR3_

#include "ito_plasma_physics_lfa.H"
#include "CD_LookupTable.H"

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace ItoPlasma{

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_air3 : public ito_plasma_physics_lfa {
    public:

      class electron;
      class positive;
      class negative;
      class PhotonZ;
      
      ito_plasma_air3();
      ~ito_plasma_air3();

      Real computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const override;

      Real computeAlpha(const RealVect a_E) const override;

      Real compute_eta(const RealVect a_E) const;

      void update_reaction_rates_lfa(const RealVect a_E, const Real a_dx, const Real a_kappa) const override;

      Vector<Real> compute_ito_mobilities_lfa(const Real a_time, const RealVect a_pos, const RealVect a_E) const override;

      Vector<Real> compute_ito_diffusion_lfa(const Real         a_time,
					     const RealVect     a_pos,
					     const RealVect     a_E,
					     const Vector<Real> a_cdr_densities) const override;


    protected:

      // for dt
      Real m_deltaX;

      // Stuff for initial conditions
      Real m_blob_radius;
      Real m_particle_weight;
      int m_num_particles;
      RealVect m_blob_center;

      // Indices for species
      int m_electron_idx;
      int m_positive_idx;
      int m_negative_idx;
      int m_PhotonZ_idx;

      // Mobile ions or not
      bool m_isMobile_ions;

      // Ion mobility
      Real m_ion_mu;
      Real m_ion_D;

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

      // Zheleznyak photoionization rates
      Real excitation_rates(const Real a_E) const;
      Real sergey_factor(const Real a_O2frac) const;
      Real photo_rate(const Real a_E) const;

      // Electron velocity
      RealVect compute_electron_velocity(const RealVect a_E) const;

      // Read file into table
      void read_tables();
    };

    
    class ito_plasma_air3::electron : public ItoSpecies {
    public:
      electron();
      ~electron();
    };

    class ito_plasma_air3::positive : public ItoSpecies {
    public:
      positive();
      ~positive();
    };

    class ito_plasma_air3::negative : public ItoSpecies {
    public:
      negative();
      ~negative();
    };

    class ito_plasma_air3::PhotonZ : public RtSpecies {
    public:
      PhotonZ();
      ~PhotonZ();

      Real getKappa(const RealVect a_pos) const override;
      
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
#include "CD_NamespaceFooter.H"

#endif
