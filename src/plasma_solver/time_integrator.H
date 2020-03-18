/*!
  @file   time_integrator.H
  @brief  Declaration of an abstract class for integrating the equations of motion
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _TIME_INTEGRATOR_
#define _TIME_INTEGRATOR_

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "amr_mesh.H"
#include "mfis.H"
#include "cdr_layout.H"
#include "rte_layout.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

#include <RealVect.H>

/*!
  @brief Abstract class for evolving plasma equations in time
*/
class time_integrator {
public:

  time_integrator(){}

  virtual ~time_integrator(){};

  virtual void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const = 0;
  //  virtual void write_plot_data() const = 0;
};

#endif
