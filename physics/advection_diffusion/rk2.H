/*!
  @file   rk2.H
  @brief  time_stepper class for only solving the Advection_Diffusion equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _ADVECTION_DIFFUSION_RK2_
#define _ADVECTION_DIFFUSION_RK2_

#include "advection_diffusion_stepper.H"

namespace physics {
  namespace advection_diffusion {

    /*!
      @brief Class that just solves the advection_diffusion equation
    */
    class rk2 : public advection_diffusion_stepper {
    public:

      rk2();
      rk2(RefCountedPtr<cdr_solver>& a_solver);
      ~rk2();

      // Pure routines for derived classes
      Real advance(const Real a_dt) final;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) final;
      void setup_solvers() final;

    protected:

      EBAMRCellData m_tmp;
      EBAMRCellData m_k1;
      EBAMRCellData m_k2;
    };
  }
}
#endif
