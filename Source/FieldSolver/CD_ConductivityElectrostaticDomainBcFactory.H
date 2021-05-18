/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ConductivityElectrostaticDomainBcFactory.H
  @brief  Contains declarations and documentation for a factory class for ConductivityElectrostaticDomainBc
  @author Robert Marskar
*/

#ifndef CD_ConductivityElectrostaticDomainBcFactory_H
#define CD_ConductivityElectrostaticDomainBcFactory_H

// Chombo includes
#include <BaseDomainBC.H>
#include <ConductivityBaseDomainBC.H>

// Our includes
#include <CD_ElectrostaticDomainBc.H>
#include <CD_ConductivityElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which interfaces a boundary condition wrapper to ConductivityBaseDomainBC.
  @details This class provides the functionality for setting Dirichlet and Neumann boundary conditions individually on each domain side, with different function. The
  boundary conditions are supplied by std::function which captures FieldSolver::m_time, so they should automatically use the correct time for the BC evaluation. 
*/
class ConductivityElectrostaticDomainBcFactory : public BaseDomainBCFactory {
public:

  /*!
    @brief Weak construction is disallowed because it muddies the interface. 
  */
  ConductivityElectrostaticDomainBcFactory() = delete;

  /*!
    @brief Full constructor which instantiates the factory. 
    @param[in] a_domainBc Domain boundary conditions.
    @param[in] a_probLo   Lower-left corner of the computational domain. 
  */
  ConductivityElectrostaticDomainBcFactory(const ElectrostaticDomainBc& a_domainBc, const RealVect a_probLo);

  /*!
    @brief Destructor. This does nothing.
  */
  ~ConductivityElectrostaticDomainBcFactory();

  /*!
    @brief Factory method. 
    @details This returns an object of the type ConductivityElectrostaticDomainBcFactory. The boundary condition functions and types are passed onto the object.
  */
  ConductivityElectrostaticDomainBc* create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx) override final;
  
protected:
  
  /*! 
    @brief Domain boundary conditions.
  */
  ElectrostaticDomainBc m_domainBc;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;
};

#include <CD_NamespaceFooter.H>

#endif
