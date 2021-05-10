/*!
  @file    dirichlet_func.H
  @brief   Declaration of a potential function class 
  @details This is the class that field_solver_multigrid uses to pass boundary conditions to AMRMultiGrid
  @author  Robert Marskar
  @date    June 2018
*/

#ifndef _DIRICHLET_FUNC_
#define _DIRICHLET_FUNC_

#include <functional>
#include <BaseBCFuncEval.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class for passing time and space dependent boundary conditions into the GMG boundary condition classes. 
*/
class dirichlet_func : public BaseBCFuncEval {
public:

  dirichlet_func() = delete;
  

  /*!
    @brief blablabla, not gonna use this anyways
  */
  dirichlet_func(std::function<Real(const Real a_time)>& a_Potential){
    m_Potential = a_Potential;
  }

  /*!
    @brief Destructor
  */
  virtual ~dirichlet_func();

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
  Real (*m_potential)(const Real a_time);

  /*!
    @brief Spatial distribution
  */
  Real (*m_func)(const RealVect a_pos);

  std::function<Real(const Real)> m_Potential;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Physical origin (the value function is called with computational coordinates)
  */
  RealVect m_origin;
};
#include "CD_NamespaceFooter.H"

#endif
