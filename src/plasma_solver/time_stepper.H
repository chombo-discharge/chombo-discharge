/*!
  @file time_stepper.H
  @brief Declaration of a time_stepper, an abstract class for evolving in time across an amr mesh
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _TIME_STEPPER_
#define _TIME_STEPPER_

#include "RealVect.H"

/*!
  @brief Abstract class for evolving in time across an amr mesh
*/
class time_stepper {
public:

  /*!
    @brief Constructor
  */
  time_stepper();

  /*!
    @brief Destructor
  */
  ~time_stepper();

  /*!
    @brief Advance method, advances equations.
    @param[in] a_dt Time step
    @returns Actual time step used (can be different from a_dt for adaptive methods)
  */
  virtual Real advance(const Real a_dt) = 0;

  /*!
    @brief Query about necessary number of ghost cells. This collects that from the individual solvers. 
  */
  virtual int query_ghost(){
    return 2;
  }

  /*!
    @brief Instantiate solvers
  */
  virtual void instantiate_solvers(){

  }
  
protected:

};

#endif
