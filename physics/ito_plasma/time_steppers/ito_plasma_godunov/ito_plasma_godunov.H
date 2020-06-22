/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_GODUNOV_
#define _ITO_PLASMA_GODUNOV_

#include "ito_plasma_stepper.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_godunov : public ito_plasma_stepper {
    public:

      // Ctor and dtor
      ito_plasma_godunov();
      ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_godunov();

      // Advance routines
      Real advance(const Real a_dt) override;

      // Setup routines
      void parse_options() override;
      void allocate_internals() override;

    protected:

      int m_merge_interval; // How often to merge superparticles

      void advect_particles(const Real a_dt);
      void diffuse_particles(const Real a_dt);
      void intersect_particles(const Real a_dt);
      void rewind_particles();
    };
  }
}
#endif
