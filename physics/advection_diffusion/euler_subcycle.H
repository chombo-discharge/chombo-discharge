/*!
  @file   euler_subcycle.H
  @brief  time_stepper class for only solving the Advection_Diffusion equation (with surface charge)
  @author Robert Marskares
  @data   March 2020
*/

#ifndef _ADVECTION_DIFFUSION_EULER_SUBCYCLE_
#define _ADVECTION_DIFFUSION_EULER_SUBCYCLE_

#include "advection_diffusion_stepper.H"

namespace physics {
  namespace advection_diffusion {

    /*!
      @brief Class that just solves the advection_diffusion equation
    */
    class euler_subcycle : public advection_diffusion_stepper {
    public:

      euler_subcycle();
      euler_subcycle(RefCountedPtr<cdr_solver>& a_solver);
      ~euler_subcycle();

      void compute_dt(Real& a_dt, time_code::which_code& a_timecode);
      Real advance(const Real a_dt) final;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) final;
      void setup_solvers() final;

    protected:

      EBAMRCellData m_k1;
      EBAMRCellData m_k2;
    };
  }
}
#endif
