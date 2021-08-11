/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzNeumannEBBCFactory.H
  @brief  Declaration of a Neumann boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzNeumannEBBCFactory_H
#define CD_MFHelmholtzNeumannEBBCFactory_H

// Our includes
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Neumann boundary conditions on the EB in a multifluid context. 
  @note This class reuses EBHelmholtzNeumannEBBCFactory, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzNeumannEBBCFactory : public MFHelmholtzEBBCFactory {
public:

  /*!
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzNeumannEBBCFactory();

  /*!
    @brief Full constructor which sets order and value on EB. 
    @param[in] a_value  Constant value on EB
    @note The order is also the attempted radius of the stencil. 
  */
  MFHelmholtzNeumannEBBCFactory(const Real a_value);

  /*!
    @brief Full constructor which sets order and spatially varying value on EBs
    @param[in] a_order Approximation order
    @param[in] a_value Spatially varying value on EBs
  */
  MFHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzNeumannEBBCFactory();

  /*!
    @brief Set value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setDphiDn(const Real a_DphiDn);

  /*!
    @brief Set spatially dependent value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setBxDphiDn(const Real a_BxDphiDn);

  /*!
    @brief Set the spatially dependent value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

  /*!
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const override;

protected:

  /*!
    @brief Multiply by b-coefficient or not
  */
  bool m_multByBco;

  /*!
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Value on EB
  */
  Real m_constantDphiDn;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif
