/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzElectrostaticDomainBCFactory.H
  @brief  Declaration of a factory class for making electrostatic domain boundary conditions with MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzElectrostaticDomainBCFactory_H
#define CD_MFHelmholtzElectrostaticDomainBCFactory_H

// Our includes
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for passing electrostatic boundary conditions into MFHelmholtzOp.
  @details This class implements MFHelmholtzDomainBCFactory and must thus return an EBHelmholtzDomainBC object for encapsulating domain boundary conditions
  in EBHelmholtzOp. This is simply a multifluid wrapper. 
*/
class MFHelmholtzElectrostaticDomainBCFactory : public MFHelmholtzDomainBCFactory
{
public:
  /*!
    @brief Weak construction is bad. Use the full constructor. 
  */
  MFHelmholtzElectrostaticDomainBCFactory() = delete;

  /*!
    @brief Full constructor
    @param[in] a_electrostaticBCs Electrostatic boundary conditions on domain sides
  */
  MFHelmholtzElectrostaticDomainBCFactory(const ElectrostaticDomainBc& a_electrostaticBCs);

  /*!
    @brief Disallowed copy constructor
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticDomainBCFactory(MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticDomainBCFactory(MFHelmholtzElectrostaticDomainBCFactory&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticDomainBCFactory&
  operator=(MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move assignment
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticDomainBCFactory&
  operator=(MFHelmholtzElectrostaticDomainBCFactory&& a_other) = delete;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~MFHelmholtzElectrostaticDomainBCFactory();

  /*!
    @brief Factory method
    @param[in] a_iphase Phase
  */
  RefCountedPtr<EBHelmholtzDomainBC>
  create(const int a_iphase) const override final;

protected:
  /*!
    @brief Domain boundary conditions
  */
  ElectrostaticDomainBc m_electrostaticBCs;
};

#include <CD_NamespaceFooter.H>

#endif
