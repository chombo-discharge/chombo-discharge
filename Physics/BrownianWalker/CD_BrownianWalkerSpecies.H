/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_BrownianWalkerSpecies.H
  @brief  Declaration of a brownian walker species for the ItoSolver
  @author Robert Marskar
*/

#ifndef CD_BrownianWalkerSpecies_H
#define CD_BrownianWalkerSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    class BrownianWalkerSpecies : public ItoSpecies {
    public:
      BrownianWalkerSpecies();
      ~BrownianWalkerSpecies();
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
      RealVect randomGaussian();
      RealVect randomDirection();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
