/*!
  @file   splitstep_forward_euler.H
  @brief  Declaration of a split step forward Euler method
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _SPLITSTEP_FORWARD_EULER_
#define _SPLITSTEP_FORWARD_EULER_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
*/
class splitstep_forward_euler : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;


};

#endif
