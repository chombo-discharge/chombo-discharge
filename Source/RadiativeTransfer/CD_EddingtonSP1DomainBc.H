/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1DomainBC.H
  @brief  Simple data holder for holding domain boundary condition types for EddingtonSP1
  @author Robert Marskar
*/

#ifndef CD_EddingtonSP1DomainBc_H
#define CD_EddingtonSP1DomainBc_H

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
  @details This class is simply meant to hold boundary conditions for use in various FieldSolver code. It provides no functionality otherwise,
  and has no dependencies, other than standard C++14. 
*/
class EddingtonSP1DomainBc {
public:

  /*!
    @brief Boundary condition type
  */
  enum class BcType {
    Dirichlet,
    Neumann,
    Larsen,
  };

  /*!
    @brief Function which maps f(R^3,t) : R. Used for setting the associated value and boundary condition type. 
  */
  using BcFunction = std::function<Real(const RealVect a_position, const Real a_time) >;

  /*!
    @brief Alias for mapping a direction and side.
  */
  using Wall = std::pair<int, Side::LoHiSide>;
  
  /*!
    @brief Alias for mapping a bctype and function. 
  */
  using Bc = std::pair<BcType, BcFunction>;
  
  /*!
    @brief Initializing constructor. This does nothing. 
  */
  EddingtonSP1DomainBc();

  /*!
    @brief Destructor. This 
  */
  virtual ~EddingtonSP1DomainBc();

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
  Bc& getBc(const Wall a_wall); 

  /*!
    @brief Returns the boundary condition for a side and direction
    @param[in] a_wall Coordinate direction and siden
    @note If the boundary condition for the associated side and direction has not been previously set, this routine will cause a run-time error. 
  */
  const Bc& getBc(const Wall a_wall) const;

protected:

  /*!
    @brief This holds the map to the boundary condition type and side. The first entry in the map is the direction and side of the boundary condition.
    The other entry in the map is the BC type, and the associated function which returns the value in space and time of the BC. 
  */
  std::map<Wall, Bc> m_bcFunctions;

};

#include <CD_NamespaceFooter.H>

#endif
