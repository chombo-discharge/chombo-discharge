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
  @brief Base class for passing domain boundary conditions into MFHelmholtzOp. 
*/
class MFHelmholtzElectrostaticDomainBCFactory : public MFHelmholtzDomainBCFactory{
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
  */
  MFHelmholtzElectrostaticDomainBCFactory(MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move constructor
  */
  MFHelmholtzElectrostaticDomainBCFactory(MFHelmholtzElectrostaticDomainBCFactory&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment
  */
  MFHelmholtzElectrostaticDomainBCFactory& operator=(MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move assignment
  */
  MFHelmholtzElectrostaticDomainBCFactory& operator=(MFHelmholtzElectrostaticDomainBCFactory&& a_other) = delete;

  /*!
    @brief Destructor
  */
  ~MFHelmholtzElectrostaticDomainBCFactory();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  MFHelmholtzElectrostaticDomainBCFactory(const MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  MFHelmholtzElectrostaticDomainBCFactory& operator=(const MFHelmholtzElectrostaticDomainBCFactory& a_other) = delete;

  /*!
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzDomainBC> create(const int a_iphase) const override final;

protected:

  /*!
    @brief Domain boundary conditions
  */
  ElectrostaticDomainBc m_electrostaticBCs;
};

#include <CD_NamespaceFooter.H>

#endif
