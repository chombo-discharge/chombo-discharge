/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LarsenCoefficients.H
  @brief  Declaration of a class that holds coefficients for Photon outflow boundary conditions
  @author Robert Marskar
*/

#ifndef CD_LarsenCoefficients_H
#define CD_LarsenCoefficients_H

// Std includes
#include <functional>

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_RobinCoefficients.H>
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that holds appropriate coefficients for Robin implementation of the Larsen boundary condition
*/
class LarsenCoefficients : public RobinCoefficients {
public:

  /*!
    @brief Don't see why anyone would need this constructor
  */
  LarsenCoefficients() = delete;

  /*!
    @brief Constructor. Source term is set to zero. 
  */
  LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
		     const Real                      a_r1,
		     const Real                      a_r2);

    /*!
    @brief Constructor. 
  */
  LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
		     const Real                      a_r1,
		     const Real                      a_r2,
		     const std::function<Real(const RealVect a_pos, const Real a_time)>& a_bcFunction);

  /*!
    @brief Destructor
  */
  ~LarsenCoefficients();

  /*!
    @brief aco. Equal to 0.5*3*kappa*kappa*(1 + 3*r2)/(1 - 2*r1)
  */
  virtual Real aco(const RealVect a_pos) const override;

  /*!
    @brief bco. Equal to -kappa
  */
  virtual Real bco(const RealVect a_pos) const override;

  /*!
    @brief Surface source term. Returns whatever is set in setRhsFunction. 
  */
  virtual Real rhs(const RealVect a_pos) const override;

  /*!
    @brief Set function for source term. 
  */
  virtual void setRhsFunction(const std::function<Real(const RealVect a_pos, const Real a_time)>& a_bcFunction);
    
protected:

  RefCountedPtr<RtSpecies> m_RtSpecies;

  Real m_reflectionCoefficientOne;
    
  Real m_reflectionCoefficientTwo;

  std::function<Real(RealVect a_pos, Real a_time) > m_bcFunction;
};

#include <CD_NamespaceFooter.H>

#endif
