/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ConductivityEddingtonSP1DomainBcFactory.H
  @brief  Contains declarations and documentation for a factory class for ConductivityEddingtonSP1DomainBc
  @author Robert Marskar
*/

#ifndef CD_ConductivityEddingtonSP1DomainBcFactory_H
#define CD_ConductivityEddingtonSP1DomainBcFactory_H

// Chombo includes
#include <BaseDomainBC.H>
#include <ConductivityBaseDomainBC.H>

// Our includes
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_ConductivityEddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which interfaces a boundary condition wrapper to ConductivityBaseDomainBC.
  @details This class provides the functionality for setting Dirichlet and Neumann boundary conditions individually on each domain side, with different function. The
  boundary conditions are supplied by std::function which captures FieldSolver::m_time, so they should automatically use the correct time for the BC evaluation. 
*/
class ConductivityEddingtonSP1DomainBcFactory : public BaseDomainBCFactory {
public:

  /*!
    @brief Weak construction is disallowed because it muddies the interface. 
  */
  ConductivityEddingtonSP1DomainBcFactory() = delete;

  /*!
    @brief Full constructor which instantiates the factory. 
    @param[in] a_domainBc Domain boundary conditions.
    @param[in] a_probLo   Lower-left corner of the computational domain. 
  */
  ConductivityEddingtonSP1DomainBcFactory(const EddingtonSP1DomainBc& a_domainBc, const RealVect a_probLo);

  /*!
    @brief Destructor. This does nothing.
  */
  ~ConductivityEddingtonSP1DomainBcFactory();

  /*!
    @brief Factory method. 
    @details This returns an object of the type ConductivityEddingtonSP1DomainBcFactory. The boundary condition functions and types are passed onto the object.
  */
  ConductivityEddingtonSP1DomainBc* create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx) override final;
  
protected:
  
  /*! 
    @brief Domain boundary conditions.
  */
  EddingtonSP1DomainBc m_domainBc;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;
};

#include <CD_NamespaceFooter.H>

#endif
