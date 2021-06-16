/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBCFactory.H
  @brief  Declaration of a factory class for EBHelmholtzNeumannEBBCFactory
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzNeumannEBBCFactory_H
#define CD_EBHelmholtzNeumannEBBCFactory_H

// Our includes
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for passing Neumann BCs into EBHelmholtzOp
*/
class EBHelmholtzNeumannEBBCFactory : public EBHelmholtzEBBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzNeumannEBBCFactory();

  /*!
    @brief Full constructor which sets dphi/dn on the EB
  */
  EBHelmholtzNeumannEBBCFactory(const Real a_DphiDn);

  /*!
    @brief Full constructor which sets dphi/dn on the EB. 
  */
  EBHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzNeumannEBBCFactory();

  /*!
    @brief Set value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setDphiDn(const int a_DphiDn);

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
  void setBxDphiDn(const int a_BxDphiDn);

  /*!
    @brief Set the spatially dependent value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

  /*!
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzEBBC> create() override;

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
