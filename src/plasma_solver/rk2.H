/*!
  @file   rk2.H
  @brief  Declaration of a second order Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _rk2_
#define _rk2_

#include "time_stepper.H"

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
*/
class rk2 : public time_stepper {
public:

  /*!
    @brief Constructor
  */
  rk2() : time_stepper() {

  }

  virtual ~rk2(){
    
  }

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt){
    return 0.;
  }

  
 protected:
};

#endif

