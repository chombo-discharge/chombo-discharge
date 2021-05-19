/*!
  @file   larsen_coefs.H
  @brief  Declaration of a class that holds coefficients for photon outflow boundary conditions
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _LARSEN_COEFS_
#define _LARSEN_COEFS_

#include <CD_RobinCoefficients.H>
#include "rte_species.H"

#include <RefCountedPtr.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class that holds appropriate coefficients for Robin implementation of the Larsen boundary condition
*/
class larsen_coefs : public RobinCoefficients {
public:

  /*!
    @brief Constructor
  */
  larsen_coefs(const RefCountedPtr<rte_species>& a_rte_species,
	       const Real                         a_r1,
	       const Real                         a_r2);

  /*!
    @brief Destructor
  */
  ~larsen_coefs();

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
#include "CD_NamespaceFooter.H"

#endif
