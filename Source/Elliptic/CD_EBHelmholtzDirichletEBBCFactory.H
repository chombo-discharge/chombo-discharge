/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBCFactory.H
  @brief  Declaration of a factory class for EBHelmholtzDirichletEBBCFactory
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletEBBCFactory_H
#define CD_EBHelmholtzDirichletEBBCFactory_H

// Our includes
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making Dirichlet BCs in EBHelmholtzEBBCFactory
*/
class EBHelmholtzDirichletEBBCFactory : public EBHelmholtzEBBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzDirichletEBBCFactory();

  /*!
    @brief Full constructor which sets order and value on EB. 
    @note The order is also the attempted radius of the stencil. 
  */
  EBHelmholtzDirichletEBBCFactory(const int a_order, const Real a_value);

  /*!
    @brief Full constructor which sets order and value on EB
  */
  EBHelmholtzDirichletEBBCFactory(const int a_order, const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDirichletEBBCFactory();

  /*!
    @brief Set BC order
    @note The order is also the attempted radius of the stencil. 
  */
  void setOrder(const int a_order);

  /*!
    @brief Set the value on the EB. 
  */
  void setValue(const Real a_value);

  /*!
    @brief Set the value on the EB. 
  */
  void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzEBBC> create() override;

protected:

  /*!
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Stencil approximation order,
  */
  int m_order;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
