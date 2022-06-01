/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzDirichletDomainBCFactory.H
  @brief  Declaration of a factory class for making Dirichlet domain boundary condition objects for MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzDirichletDomainBCFactory_H
#define CD_MFHelmholtzDirichletDomainBCFactory_H

// Std includes
#include <functional>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for passing Dirichlet domain boundary conditions into MFHelmholtzOp. 
*/
class MFHelmholtzDirichletDomainBCFactory : public MFHelmholtzDomainBCFactory
{
public:
  /*!
    @brief Constructor. 
  */
  MFHelmholtzDirichletDomainBCFactory();

  /*!
    @brief Constructor. 
    @param[in] a_value Constant value on domain side
  */
  MFHelmholtzDirichletDomainBCFactory(const Real a_value);

  /*!
    @brief Full constructor. 
    @param[in] a_value Spatially dependent value on domain side
  */
  MFHelmholtzDirichletDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  virtual ~MFHelmholtzDirichletDomainBCFactory();

  /*!
    @brief Set the constant value on the domain side. 
    @param[in] a_value Constant value on all domain sides. 
  */
  void
  setValue(const Real a_value);

  /*!
    @brief Set the spatially dependent value on the domain side. 
    @param[in] a_value Spatially dependent value on all domain sides. 
  */
  void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Factory method
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC>
  create(const int a_iphase) const override;

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
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
