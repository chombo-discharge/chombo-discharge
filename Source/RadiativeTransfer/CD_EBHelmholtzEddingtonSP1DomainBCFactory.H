/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEddingtonSP1DomainBCFactory.H
  @brief  Declaration of a factory class for wrapping domain boundary conditions for multigrid in the EddingtonSP1 solver
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEddingtonSP1DomainBCFactory_H
#define CD_EBHelmholtzEddingtonSP1DomainBCFactory_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EddingtonSP1 domain boundary conditions into EBHelmholtzOp
*/
class EBHelmholtzEddingtonSP1DomainBCFactory : public EBHelmholtzDomainBCFactory
{
public:
  /*!
    @brief Weak construction is bad. Use the full constructor. 
  */
  EBHelmholtzEddingtonSP1DomainBCFactory() = delete;

  /*!
    @brief Full constructor
    @param[in] a_eddingtonBCs EddingtonSP1 boundary conditions on domain sides
    @param[in] a_species      Radiative transfer species
    @param[in] a_r1           "First reflection coefficient"
    @param[in] a_r2           "Second reflection coefficient"
  */
  EBHelmholtzEddingtonSP1DomainBCFactory(const EddingtonSP1DomainBc&     a_eddingtonBCs,
                                         const RefCountedPtr<RtSpecies>& a_species,
                                         const Real                      a_r1,
                                         const Real                      a_r2);

  /*!
    @brief Disallowed copy constructor
  */
  EBHelmholtzEddingtonSP1DomainBCFactory(EBHelmholtzEddingtonSP1DomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move constructor
  */
  EBHelmholtzEddingtonSP1DomainBCFactory(EBHelmholtzEddingtonSP1DomainBCFactory&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment
  */
  EBHelmholtzEddingtonSP1DomainBCFactory&
  operator=(EBHelmholtzEddingtonSP1DomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed move assignment
  */
  EBHelmholtzEddingtonSP1DomainBCFactory&
  operator=(EBHelmholtzEddingtonSP1DomainBCFactory&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzEddingtonSP1DomainBCFactory();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzEddingtonSP1DomainBCFactory(const EBHelmholtzEddingtonSP1DomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzEddingtonSP1DomainBCFactory&
  operator=(const EBHelmholtzEddingtonSP1DomainBCFactory& a_other) = delete;

  /*!
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzDomainBC>
  create() const override final;

protected:
  /*!
    @brief Domain boundary conditions
  */
  EddingtonSP1DomainBc m_eddingtonBCs;

  /*!
    @brief Radiative transfer species
  */
  RefCountedPtr<RtSpecies> m_species;

  /*!
    @brief "First reflection coefficient"
  */
  Real m_r1;

  /*!
    @brief "Second reflection coefficient"
  */
  Real m_r2;
};

#include <CD_NamespaceFooter.H>

#endif
