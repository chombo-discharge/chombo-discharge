/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RobinConductivityEbBcFactory.H
  @brief  Declaration of a factory class for Robin type EB boundary condition
  @author Robert Marskar
*/

#ifndef CD_RobinConductivityEbBcFactory_H
#define CD_RobinConductivityEbBcFactory_H

// Chombo includes
#include <BaseEBBC.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_RobinCoefficients.H>
#include <CD_IrregStencil.H>
#include <CD_RobinConductivityEbBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for generating Robin boundary conditions \f$a\phi + b\partial_n\phi = g\f$, where \f$n\f$ is the EB unit normal. 
*/
class RobinConductivityEbBcFactory : public BaseEBBCFactory {
public:

  /*!
    @brief Constructor
  */
  RobinConductivityEbBcFactory(const RealVect a_origin);

  /*!
    @brief Destructor
  */
  ~RobinConductivityEbBcFactory();

  /*!
    @brief Set bc order
  */
  virtual void setStencilType(const IrregStencil::StencilType a_type);

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
  virtual void setCoefficients(const EBAMRIVData& a_aco,
			       const EBAMRIVData& a_bco,
			       const EBAMRIVData& a_rhs);

  /*!
    @brief Create Robin type boundary condition
  */
  RobinConductivityEbBc* create(const ProblemDomain& a_domain,
				const EBISLayout&    a_ebisl,
				const RealVect&      a_dx,
				const IntVect*       a_ghost_phi = 0,
				const IntVect*       a_ghost_rhs = 0);

protected:

  IrregStencil::StencilType m_type;

  RealVect m_origin;

  RefCountedPtr<RobinCoefficients> m_robinCoefficients;

  Real m_aco;

  Real m_bco;
  
  Real m_rhs;

  EBAMRIVData m_acodata;

  EBAMRIVData m_bcodata;
  
  EBAMRIVData m_rhsdata;

  bool m_const_coeff;

  bool m_data_coeff;

  bool m_func_coeff;
};

#include <CD_NamespaceFooter.H>

#endif
