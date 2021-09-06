/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBCFactory.H
  @brief  Declaration of a factory class for making Neumann domain boundary condition objects for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzNeumannDomainBCFactory_H
#define CD_EBHelmholtzNeumannDomainBCFactory_H

// Std includes
#include <functional>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for passing Neumann boundary domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzNeumannDomainBCFactory : public EBHelmholtzDomainBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzNeumannDomainBCFactory();

  /*!
    @brief Constructor. 
    @param[in] a_DphiDn Constant value of dphi/dn on all domain sides. 
    @note A positive values gives a flux INTO the domain. 
  */
  EBHelmholtzNeumannDomainBCFactory(const Real a_DphiDn);

  /*!
    @brief Full constructor. 
    @param[in] a_DphiDn Spatially dependent value of dphi/dn on all domain sides. 
    @note A positive values gives a flux INTO the domain. 
  */
  EBHelmholtzNeumannDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzNeumannDomainBCFactory();

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
    @brief Factory method. 
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC> create() const override;

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
    @brief Constant values on domain sides. 
  */
  Real m_constantDphiDn;

  /*!
    @brief Spatially dependent value on the domain sides. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionDphiDn;
  
};

#include <CD_NamespaceFooter.H>

#endif
