/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RobinConductivityDomainBcFactory.H
  @brief  Declaration of a factory class for Robin boundary conditions on domain faces
  @author Robert Marskar
*/

#ifndef CD_RobinConductivityDomainBcFactory_H
#define CD_RobinConductivityDomainBcFactory_H

// Chombo includes
#include <BaseDomainBC.H>

// Our includes
#include <EBAMRData.H>
#include <CD_RobinCoefficients.H>
#include <CD_RobinConductivityDomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for generating Robin boundary condition on domain faces
*/
class RobinConductivityDomainBcFactory : public BaseDomainBCFactory {
public:

  /*!
    @brief Constructor
  */
  RobinConductivityDomainBcFactory();

  /*!
    @brief Destructor
  */
  ~RobinConductivityDomainBcFactory();

  /*!
    @brief Function based coefficients and right hand sides
  */
  virtual void setCoefficients(RefCountedPtr<RobinCoefficients> a_robinco);

  /*!
    @brief Set constant coefficients
  */
  virtual void setCoefficients(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Set data-based coefficients
  */
  virtual void setCoefficients(const EBAMRFluxData& a_aco, const EBAMRFluxData& a_bco, const EBAMRFluxData& a_rhs);

  /*!
    @brief Create robin boundary condition object
  */
  virtual BaseDomainBC* create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx);
  
protected:

  RefCountedPtr<RobinCoefficients> m_robinco;

  Real m_aco;

  Real m_bco;
  
  Real m_rhs;

  EBAMRFluxData m_acodata;

  EBAMRFluxData m_bcodata;
  
  EBAMRFluxData m_rhsdata;

  bool m_const_coeff;

  bool m_data_coeff;

  bool m_func_coeff;
};

#include <CD_NamespaceFooter.H>

#endif
