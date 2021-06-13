/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBCFactory.H
  @brief  Declaration of a factory class for making boundary condition objects for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletDomainBCFactory_H
#define CD_EBHelmholtzDirichletDomainBCFactory_H

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
class EBHelmholtzDirichletDomainBCFactory : public EBHelmholtzDomainBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzDirichletDomainBCFactory();

  /*!
    @brief Constructor. 
    @param[in] a_value Constant value on domain side
  */
  EBHelmholtzDirichletDomainBCFactory(const Real a_value);

  /*!
    @brief Constructor. 
    @param[in] a_value Spatially dependent value on domain side
  */
  EBHelmholtzDirichletDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDirichletDomainBCFactory();

  /*!
    @brief Set the value on the domain side. 
  */
  void setValue(const Real a_value);

  /*!
    @brief Set the value on the domain side. 
  */
  void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Factory method
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC> create() override;

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
