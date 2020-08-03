/*!
  @file   robinconductivitydomainbcfactory.H
  @brief  Declaration of a factory class for Robin boundary conditions on domain faces
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _ROBINCONDUCTIVITYDOMAINBCFACTORY_
#define _ROBINCONDUCTIVITYDOMAINBCFACTORY_

#include "EBAMRData.H"
#include "robin_coef.H"
#include "robinconductivitydomainbc.H"

#include <BaseDomainBC.H>

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
  virtual void set_coefs(RefCountedPtr<robin_coef> a_robinco);

  /*!
    @brief Set constant coefficients
  */
  virtual void set_coefs(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Set data-based coefficients
  */
  virtual void set_coefs(const EBAMRFluxData& a_aco, const EBAMRFluxData& a_bco, const EBAMRFluxData& a_rhs);

  /*!
    @brief Create robin boundary condition object
  */
  virtual BaseDomainBC* create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx);
  
protected:

  RefCountedPtr<robin_coef> m_robinco;

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

#endif
