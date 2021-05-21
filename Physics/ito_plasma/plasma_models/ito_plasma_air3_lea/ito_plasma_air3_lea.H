/*!
  @file   ito_plasma_air3_lea.H
  @brief  Simple class for 3-species air plasma
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_AIR3_LEA_
#define _ITO_PLASMA_AIR3_LEA_

#include "ito_plasma_physics_lea.H"
#include "CD_LookupTable.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace ito_plasma{

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_air3_lea : public ito_plasma_physics_lea {
    public:

      class electron;
      class positive;
      class negative;
      class PhotonZ;
      
      ito_plasma_air3_lea();
      ~ito_plasma_air3_lea();

      Real computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const override;      

      Real compute_alpha(const RealVect a_E) const override;

      Real compute_eta(const RealVect a_E) const;

      void update_reaction_rates_lea(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const override;

    protected:

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


      // Pressure, neutral number density etc
      Real m_N;      
      Real m_p;      
      Real m_pq;     
      Real m_T;      
      Real m_O2frac; 
      Real m_N2frac; 
      Real m_pO2;    
      Real m_photoi_factor;

      // Zheleznyak photoionization rates
      Real excitation_rates(const Real a_E) const;
      Real sergey_factor(const Real a_O2frac) const;
      Real photo_rate(const Real a_E) const;

      // Read file into table
      void read_tables();
    };

    
    class ito_plasma_air3_lea::electron : public ItoSpecies {
    public:
      electron() = default;
      electron(const LookupTable& a_mobility, const LookupTable& a_diffusion);
      ~electron();

      Real mobility(const Real a_energy) const override;
      Real diffusion(const Real a_energy) const override;

    protected:
      LookupTable m_mobility;
      LookupTable m_diffusion;
    };

    class ito_plasma_air3_lea::positive : public ItoSpecies {
    public:
      positive();
      ~positive();
    };

    class ito_plasma_air3_lea::negative : public ItoSpecies {
    public:
      negative();
      ~negative();
    };

    class ito_plasma_air3_lea::PhotonZ : public RtSpecies {
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
