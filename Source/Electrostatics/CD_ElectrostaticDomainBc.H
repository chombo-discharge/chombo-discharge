/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ElectrostaticDomainBc.H
  @brief  Declaration of a class for holding function-based boundary conditions on domain edges. 
  @author Robert Marskar
*/

#ifndef CD_ElectrostaticDomainBc_H
#define CD_ElectrostaticDomainBc_H

// Std include
#include <functional>
#include <map>

// Chombo includes
#include <RealVect.H>
#include <LoHiSide.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which maps boundary condition types to a side and direction. 
  @details This class is simply meant to hold functions that will be used as boundary conditions for use in various FieldSolver code. The user
  will construct this object, after which he can associate a polymorphic function with a domain side and a boundary condition type. Currently,
  only Dirichlet and Neumann are supported. 
*/
class ElectrostaticDomainBc
{
public:
  /*!
    @brief Supported boundary condition types.
  */
  enum class BcType
  {
    Dirichlet,
    Neumann
  };

  /*!
    @brief Function which maps f(R^3,t) : R. Used for setting the associated value and boundary condition type. 
    @param[in] a_position Physical coordinates
    @param[in] a_time     Time
  */
  using BcFunction = std::function<Real(const RealVect a_position, const Real a_time)>;

  /*!
    @brief Alias for mapping a direction and side.
  */
  using DomainSide = std::pair<int, Side::LoHiSide>;

  /*!
    @brief Alias for mapping a bctype and function. 
  */
  using Bc = std::pair<BcType, BcFunction>;

  /*!
    @brief Initializing constructor. This populates the BC functions with homogeneous Neumann boundary conditions. 
  */
  ElectrostaticDomainBc();

  /*!
    @brief Copy constructor
    @param[in] a_other Other
  */
  ElectrostaticDomainBc(const ElectrostaticDomainBc& a_other) = default;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other
  */
  ElectrostaticDomainBc(const ElectrostaticDomainBc&& a_other) = delete;

  /*!
    @brief Copy assignment
    @param[in] a_other Other
  */
  ElectrostaticDomainBc&
  operator=(const ElectrostaticDomainBc& a_other) = default;

  /*!
    @brief Disallowed move assignment
    @param[in] a_other Other
  */
  ElectrostaticDomainBc&
  operator=(const ElectrostaticDomainBc&& a_other) = delete;

  /*!
    @brief Destructor. This 
  */
  virtual ~ElectrostaticDomainBc();

  /*!
    @brief Set the boundary condition on a domain side. You must set both the type, and the associated function. 
    @param[in] a_domainSide Domain side
    @param[in] a_bc         BC type and function
    @note If the boundary condition had been set previously, this code will overwrite the previous setting.
  */
  void
  setBc(const DomainSide a_wall, const Bc a_bc);

  /*!
    @brief Returns the boundary condition for a domain side
    @param[in] a_domainSide Domain side
    @note If the boundary condition for the associated side and direction has not been previously set, this routine will cause a run-time error. 
  */
  Bc
  getBc(const DomainSide a_domainSide) const;

protected:
  /*!
    @brief This holds the map to the boundary condition type and side. The first entry in the map is the direction and side of the boundary condition.
    The other entry in the map is the BC type, and the associated function which returns the value in space and time of the BC. 
  */
  std::map<DomainSide, Bc> m_bcFunctions;
};

#include <CD_NamespaceFooter.H>

#endif
