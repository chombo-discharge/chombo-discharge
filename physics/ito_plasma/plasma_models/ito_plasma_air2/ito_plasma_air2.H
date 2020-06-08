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

    protected:

      // Stuff for initial conditions
      Real     m_blob_radius;
      int      m_num_particles;
      RealVect m_blob_center;

      int m_seed;
      std::mt19937_64 m_rng;
      std::uniform_real_distribution<Real> m_udist11;
      std::normal_distribution<Real> m_gauss;

      void draw_initial_particles();
      RealVect random_gaussian();
      RealVect random_direction();

      int m_electron_idx;
      int m_positive_idx;

      
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

    
  }
}
#endif
