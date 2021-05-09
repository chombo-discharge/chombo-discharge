/*!
  @file   CD_ElectrostaticDomainBc.H
  @brief  Simple data holder for holding domain boundary condition types for field_solver
  @author Robert Marskar
  @date   May 2021
*/

#ifndef _CD_ElectrostaticDomainBc_H_
#define _CD_ElectrostaticDomainBc_H_

#include <functional>
#include <map>

#include "RealVect.H"
#include "LoHiSide.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Class which maps boundary condition types to a side and direction. 
  @details This class is simply meant to hold boundary conditions for use in various field_solver code. It provides no functionality otherwise,
  and has no dependencies, other than standard C++14. 
*/
class ElectrostaticDomainBc {
public:

  /*!
    @brief Boundary condition type
  */
  enum class BcType {
    Dirichlet,
    Neumann
  };

  /*!
    @brief Function which maps f(R^3,t) : R. Used for setting the associated value and boundary condition type. 
  */
  using BcFunction = std::function<Real(const RealVect a_position, const Real a_time) >;

  /*!
    @brief Shortcut for mapping a direction and side.
  */
  using Wall = std::pair<int, Side::LoHiSide>;
  
  /*!
    @brief Shortcut for mapping a bctype and function. 
  */
  using Bc = std::pair<BcType, BcFunction>;
  
  /*!
    @brief Initializing constructor. This does nothing. 
  */
  ElectrostaticDomainBc();

  /*!
    @brief Destructor. 
  */
  virtual ~ElectrostaticDomainBc();

  /*!
    @brief Set the boundary condition on a domain side. You must set both the type, and the associated function. 
    @param[in] a_wall Coordinate direction and side
    @param[in] a_bc   BC type and function
    @note If the boundary condition had been set previously, this code will overwrite the previous setting.
  */  
  void setBc(const Wall a_wall, const Bc a_func);

  /*!
    @brief Returns the boundary condition for a side and direction
    @param[in] a_wall Coordinate direction and siden
    @note If the boundary condition for the associated side and direction has not been previously set, this routine will cause a run-time error. 
  */
  Bc getBc(const Wall a_wall) const;

protected:

  /*!
    @brief This holds the map to the boundary condition type and side. The first entry in the map is the direction and side of the boundary condition.
    The other entry in the map is the BC type, and the associated function which returns the value in space and time of the BC. 
  */
  std::map<Wall, Bc> m_bcFunctions;

};

#include "CD_NamespaceFooter.H"

#endif
