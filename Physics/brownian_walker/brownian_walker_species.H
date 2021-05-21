/*!
  @file   brownian_walker_species.H
  @brief  Declaration of a brownian walker species for the ito_solver
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _BROWNIAN_WALKER_SPECIES_
#define _BROWNIAN_WALKER_SPECIES_

#include <CD_ItoSpecies.H>

#include <random>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace brownian_walker {

    class brownian_walker_species : public ItoSpecies {
    public:
      brownian_walker_species();
      ~brownian_walker_species();
    protected:
      Real     m_blob_radius;
      Real     m_blob_amplitude;
      int      m_num_particles;
      RealVect m_blob_center;

      int m_seed;
      std::mt19937_64 m_rng;
      std::uniform_real_distribution<Real> m_udist11;
      std::normal_distribution<Real> m_gauss;

      void draw_initial_particles();
      RealVect random_gaussian();
      RealVect randomDirection();
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
