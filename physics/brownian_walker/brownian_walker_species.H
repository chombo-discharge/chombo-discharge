/*!
  @file   brownian_walker_species.H
  @brief  Declaration of a brownian walker species for the ito_solver
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _BROWNIAN_WALKER_SPECIES_
#define _BROWNIAN_WALKER_SPECIES_

#include "ito_species.H"

namespace physics {
  namespace brownian_walker {

    class brownian_walker_species : public ito_species {
    public:
      brownian_walker_species();
      ~brownian_walker_species();
    protected:
      Real     m_blob_radius;
      Real     m_blob_amplitude;
      int      m_num_particles;
      RealVect m_blob_center;
    };
  }
}
#endif
