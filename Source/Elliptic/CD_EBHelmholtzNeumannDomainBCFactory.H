/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBCFactory.H
  @brief  Declaration of a factory class for making boundary condition objects for EBHelmholtzOp
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
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzNeumannDomainBCFactory : public EBHelmholtzDomainBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzNeumannDomainBCFactory();

  /*!
    @brief Constructor. 
    @param[in] a_value Constant value on domain side
  */
  EBHelmholtzNeumannDomainBCFactory(const Real a_DphiDn);

  /*!
    @brief Constructor. 
    @param[in] a_value Spatially dependent value on domain side
  */
  EBHelmholtzNeumannDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzNeumannDomainBCFactory();

  /*!
    @brief Set the value of dphi/dn on the EB. 
  */
  void setDphiDn(const int a_DphiDn);

  /*!
    @brief Set the value of dphi/dn on the EB. 
  */
  void setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
  */
  void setBxDphiDn(const int a_BxDphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
  */
  void setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);
  
  /*!
    @brief Factory method
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC> create() override;

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
