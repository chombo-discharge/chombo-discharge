/*!
  @file   CD_ElectrostaticDomainBcFuncEval.H
  @brief  Contains declaration of a middle-man class for parsing domain boundary conditions into a language understandable for Chombo multigrid
  @author Robert Marskar
  @date   May 2021
*/

#ifndef _CD_ElectrostaticDomainBcFuncEval_H_
#define _CD_ElectrostaticDomainBcFuncEval_H_

#include <functional>

#include "BaseBCFuncEval.H"

#include "CD_ElectrostaticDomainBc.H"
#include "CD_NamespaceHeader.H"

/*!
  @brief Class which interfaces from more general boundary condition specification to something usable by the current Chombo multigrid solvers that field_solver_multigrid uses. 
*/
class ElectrostaticDomainBcFuncEval : public BaseBCFuncEval {
public:
  
  /*!
    @brief Deleted initializing constructor. 
  */
  ElectrostaticDomainBcFuncEval() = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_bcFunc The boundary condition function.
    @param[in] a_probLo Lower-left corner of computational domain. 
    @note The origin is sent in because EBConductivityOp calls the boundary condition function with no origin. So, this evaluation of BaseBCFuncEval rectifies that.
  */
  ElectrostaticDomainBcFuncEval(const ElectrostaticDomainBc::BcFunction a_bcFunc, const RealVect a_probLo);

  /*!
    @brief Destructor. This does nothing. 
  */
  ~ElectrostaticDomainBcFuncEval();

  /*!
    @brief Set time
    @param[in] a_time Time for when the bcFunc value is updated. 
  */
  void setTime(const Real a_time);

  /*!
    @brief Return functional evaluation in space and time of the boundary condition.
    @param[in] a_point Point in space
    @param[in] a_comp  Component. 
    @details The time-dependence is implicit - the user needs to set the time by calling setTime() first. 
    @note Calling this without having set the function first yields an error. 
  */
  Real value(const RealVect& a_point, const int& a_comp) const override final;

  /*!
    @brief Implementation of BaseBCFuncEval.
    @param[in] a_point    Position in space
    @param[in] a_comp     Component in elliptic equation.
    @param[in] a_derivDir Derivative direction.
    @return Throws an error if called since this function is not used in field_solver_multigrid code. 
  */
  Real derivative(const RealVect& a_point, const int& a_comp, const int& a_derivDir) const override final;

protected:

  /*!
    @brief Time for the domain boundary function
  */
  Real m_time;
  
  /*!
    @brief Lower left corner of computational boundary. 
  */
  RealVect m_probLo;

  /*!
    @brief Space- and time-dependent boundary condition function. 
    @note This is just an alias. 
  */
  ElectrostaticDomainBc::BcFunction m_bcFunc;
};

#include "CD_NamespaceFooter.H"

#endif
