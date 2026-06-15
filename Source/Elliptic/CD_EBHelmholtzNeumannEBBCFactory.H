/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_EBHelmholtzNeumannEBBCFactory.H
  @brief  Declaration of a factory class for making Neumann BCs in an EBHelmholtzOp context
  @author Robert Marskar
*/

#ifndef CD_EBHELMHOLTZNEUMANNEBBCFACTORY_H
#define CD_EBHELMHOLTZNEUMANNEBBCFACTORY_H

// Our includes
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Factory class which generates EBHelmholtzNeumannEBBC. 
*/
class EBHelmholtzNeumannEBBCFactory : public EBHelmholtzEBBCFactory
{
public:
  /**
    @brief Constructor. 
  */
  EBHelmholtzNeumannEBBCFactory();

  /**
    @brief Full constructor which sets dphi/dn on the EB
    @param[in] a_DphiDn Description
  */
  EBHelmholtzNeumannEBBCFactory(const Real a_DphiDn);

  /**
    @brief Full constructor which sets dphi/dn on the EB. 
    @param[in] a_DphiDn Description
  */
  EBHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /**
    @brief Destructor
  */
  virtual ~EBHelmholtzNeumannEBBCFactory();

  /**
    @brief Set value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const Real a_DphiDn);

  /**
    @brief Set spatially dependent value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /**
    @brief Set the value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const Real a_BxDphiDn);

  /**
    @brief Set the spatially dependent value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

  /**
    @brief Factory method
    @return Return value
  */
  RefCountedPtr<EBHelmholtzEBBC>
  create() override;

protected:
  /**
    @brief Multiply by b-coefficient or not
  */
  bool m_multByBco;

  /**
    @brief Use constant for BC
  */
  bool m_useConstant;

  /**
    @brief Use function for BC value
  */
  bool m_useFunction;

  /**
    @brief Value on EB
  */
  Real m_constantDphiDn;

  /**
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif
