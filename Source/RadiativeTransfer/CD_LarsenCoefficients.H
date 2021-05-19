/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_LarsenCoefficients.H
  @brief  Declaration of a class that holds coefficients for Photon outflow boundary conditions
  @author Robert Marskar
*/

#ifndef CD_LarsenCoefficients_H
#define CD_LarsenCoefficients_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_RobinCoefficients.H>
#include <rte_species.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that holds appropriate coefficients for Robin implementation of the Larsen boundary condition
*/
class LarsenCoefficients : public RobinCoefficients {
public:

  /*!
    @brief Constructor
  */
  LarsenCoefficients(const RefCountedPtr<rte_species>& a_rte_species,
		     const Real                         a_r1,
		     const Real                         a_r2);

  /*!
    @brief Destructor
  */
  ~LarsenCoefficients();

  /*!
    @brief aco. Equal to 0.5*3*kappa*kappa*(1 + 3*r2)/(1 - 2*r1)
  */
  virtual Real aco(const RealVect a_pos) const;

  /*!
    @brief bco. Equal to -kappa
  */
  virtual Real bco(const RealVect a_pos) const;

  /*!
    @brief Surface source term. Return zero.
  */
  virtual Real rhs(const RealVect a_pos) const;
    
protected:

  RefCountedPtr<rte_species> m_rte_species;

  Real m_reflectionCoefficientOne;
    
  Real m_reflectionCoefficientTwo;
};

#include <CD_NamespaceFooter.H>

#endif
