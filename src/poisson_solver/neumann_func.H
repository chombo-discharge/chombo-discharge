/*!
  @file    neumann_func.H
  @brief   Declaration of a potential function class 
  @details This is the class that poisson_multifluid_gmg uses to pass boundary conditions to AMRMultiGrid
  @author  Robert Marskar
  @date    June 2018
*/

#ifndef _NEUMANN_FUNC_
#define _NEUMANN_FUNC_

#include <BaseBCFuncEval.H>

/*!
  @brief Class for passing time and space dependent boundary conditions into the GMG boundary condition classes. 
*/
class neumann_func : public BaseBCFuncEval {
public:

  /*!
    @brief New constructor
  */
  neumann_func(Real (*a_value)(const Real a_time), Real (*a_func)(const RealVect a_pos), const RealVect a_origin);

  /*!
    @brief Destructor
  */
  virtual ~neumann_func();

  /*!
    @brief Update the time
  */
  virtual void set_time(const Real a_time);

  /*!
    @brief Value function. Called by BC classes in AMRMultiGrid
  */
  virtual Real value(const RealVect& a_point, const int& a_comp) const;

  /*!
    @brief Should never be called.
  */
  virtual Real derivative(const RealVect& a_point, const int& a_comp, const int& a_dir) const;

protected:

  /*!
    @brief Potential as function of time
  */
  Real (*m_value)(const Real a_time);

  /*!
    @brief Spatial distribution
  */
  Real (*m_func)(const RealVect a_pos);

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Physical origin (the value function is called with computational coordinates)
  */
  RealVect m_origin;
};

#endif
