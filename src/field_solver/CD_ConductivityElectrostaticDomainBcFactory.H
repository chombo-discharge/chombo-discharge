/*!
  @file   CD_ConductivityElectrostaticDomainBcFactory.H
  @brief  Contains declarations and documentation for a factory class for ConductivityElectrostaticDomainBc
  @author Robert Marskar
  @date   May 2021
*/

#ifndef _CD_ConductivityElectrostaticDomainBcFactory_H_
#define _CD_ConductivityElectrostaticDomainBcFactory_H_

#include <memory>

#include <BaseDomainBC.H>
#include <ConductivityBaseDomainBC.H>

#include "CD_ElectrostaticDomainBc.H"
#include "CD_ElectrostaticDomainBcFuncEval.H"
#include "CD_ConductivityElectrostaticDomainBc.H"
#include "CD_NamespaceHeader.H"

/*!
  @brief Class which interfaces a boundary condition wrapper to ConductivityBaseDomainBC.
  @details This class provides the functionality for setting Dirichlet and Neumann boundary conditions individually on each domain side, with different function. The
  boundary conditions are supplied by reference, so that the user can update the time for which the boundary condition is evaluated (through ElectrostaticDomainBcFuncEval).
*/
class ConductivityElectrostaticDomainBcFactory : public BaseDomainBCFactory {
public:

  /*!
    @brief Alias for easily associating a function (a BaseBCFuncEval) with a domain side.
  */
  using WallBcFuncs = std::map<ElectrostaticDomainBc::Wall, std::shared_ptr<ElectrostaticDomainBcFuncEval> >;

  /*!
    @brief Alias for easily associating a boundary condition type (an ElectrostaticDomainBC::BcType) with a domain side.
  */
  using WallBcTypes = std::map<ElectrostaticDomainBc::Wall, ElectrostaticDomainBc::BcType>;

  /*!
    @brief Constructor. This does nothing. 
  */
  ConductivityElectrostaticDomainBcFactory() = delete;

  /*!
    @brief Full constructor which instantiates the factory. 
    @param[in] a_bcFunctions Boundary condition functions on each domain side. 
    @param[in] a_bcTypes     Boundary condition types on each domain side. 
    @details The injection of the function on the domain sides is done by pointer since field_solvers may choose to use different time-dependence for his boundary conditions. 
    @param[in] a_electrostaticBCs Specification of boundary conditions. 
  */
  ConductivityElectrostaticDomainBcFactory(const WallBcFuncs& a_bcFunctions, const WallBcTypes& a_bcTypes);

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
    @brief Boundary condition functions
  */
  WallBcFuncs m_bcFunctions;

  /*!
    @brief Boundary condition types
  */
  WallBcTypes m_bcTypes;
};

#include "CD_NamespaceFooter.H"

#endif
