/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrDomainBC.H
  @brief  Declaration of domain boundary condition types for CdrSolver
  @author Robert Marskar
*/

#ifndef CD_CdrDomainBC_H
#define CD_CdrDomainBC_H

// Std includes
#include <map>
#include <functional>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for encapsulating possible domain boundary conditions for CdrSolver
  @note This class is used in CdrSolver to populate boundary conditions, letting the solver pick between various 
  scenarios, such as wall boundary conditions, data-based boundary conditions, function-based boundary conditions etc. This
  class is just a placeholder for that functionality -- it does not provide any discretization whatsoever. 
  @note This class will be run-time safe in the sense that its default constructor provides a fully functional object which (shouldn't) 
  cause run-time errors. 
*/
class CdrDomainBC {
public:
  
  /*!
    @brief Boundary condition types for the CDR solvers
  */
  enum class BcType {
    DataBased,
    Wall,
    Function,
    Outflow,
    External
  };

  /*!
    @brief Function which maps f(R^3,t) : R. Used for setting the associated value and boundary condition type. 
    @param[in] a_position Physical position
    @param[in] a_time     Time
    @return Returns flux in time and space
  */
  using FluxFunction = std::function<Real(const RealVect a_position, const Real a_time) >;

  /*!
    @brief Alias for mapping a direction and side. The first element in pair is the coordinate direction and the second is the low/high end in that direction. 
  */
  using DomainSide = std::pair<int, Side::LoHiSide>;

  /*!
    @brief Initializing constructor.
    @details This populates m_bcTypes and m_bcFunctions with BcType::Wall and flux functions which return zero (i.e., also wall). 
  */
  CdrDomainBC();

  /*!
    @brief Destructor. Does nothing. 
  */
  ~CdrDomainBC();

  /*!
    @brief Get the BC type on a particular domain side. 
    @param[in] a_domainSide Domain side
    @returns Returns the BC type on the domain side
  */
  const BcType& getBcType(const DomainSide& a_domainSide) const;

    /*!
    @brief Get the BC type on a particular domain side
    @param[in] a_domainSide Domain side
    @returns Returns the BC function on the domain side
  */
  const FluxFunction& getBcFunction(const DomainSide& a_domainSide) const;

  /*!
    @brief Set the BC type on a domain side
    @param[in] a_domainSide Domain side
    @param[in] a_bcType     Domain BC type
  */
  void setBcType(const DomainSide& a_domainSide, const BcType& a_bcType);

  /*!
    @brief Set the BC type on a domain side. 
    @param[in] a_domainSide Domain side
    @param[in] a_function   Flux function on domain side
  */
  void setBcFunction(const DomainSide& a_domainSide, const FluxFunction& a_function);

protected:

  /*!
    @brief BC types on various domain edges
  */
  std::map<DomainSide, BcType> m_bcTypes;

  /*!
    @brief BC functions on various domain edges
  */
  std::map<DomainSide, FluxFunction> m_bcFunctions;
};

#include <CD_NamespaceFooter.H>

#endif
