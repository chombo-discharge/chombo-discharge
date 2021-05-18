/*!
  @file   robinconductivitydomainbcfactory.H
  @brief  Declaration of a factory class for Robin boundary conditions on domain faces
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _ROBINCONDUCTIVITYDOMAINBCFACTORY_
#define _ROBINCONDUCTIVITYDOMAINBCFACTORY_

#include "EBAMRData.H"
#include <CD_RobinCoefficients.H>
#include "robinconductivitydomainbc.H"

#include <BaseDomainBC.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Factory class for generating Robin boundary condition on domain faces
*/
class robinconductivitydomainbcfactory : public BaseDomainBCFactory {
public:

  /*!
    @brief Constructor
  */
  robinconductivitydomainbcfactory();

  /*!
    @brief Destructor
  */
  ~robinconductivitydomainbcfactory();

  /*!
    @brief Function based coefficients and right hand sides
  */
  virtual void setCoefficientss(RefCountedPtr<RobinCoefficients> a_robinco);

  /*!
    @brief Set constant coefficients
  */
  virtual void setCoefficientss(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Set data-based coefficients
  */
  virtual void setCoefficientss(const EBAMRFluxData& a_aco, const EBAMRFluxData& a_bco, const EBAMRFluxData& a_rhs);

  /*!
    @brief Create robin boundary condition object
  */
  virtual BaseDomainBC* create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx);
  
protected:

  RefCountedPtr<RobinCoefficients> m_robinco;

  Real m_aCoefficient;

  Real m_bco;
  
  Real m_rhs;

  EBAMRFluxData m_aCoefficientdata;

  EBAMRFluxData m_bcodata;
  
  EBAMRFluxData m_rhsdata;

  bool m_const_coeff;

  bool m_data_coeff;

  bool m_func_coeff;
};
#include "CD_NamespaceFooter.H"

#endif
