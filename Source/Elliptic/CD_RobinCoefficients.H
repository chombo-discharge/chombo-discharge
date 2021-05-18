/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RobinCoefficients.H
  @brief  Declaration of a Robin coefficients for BC interfacing
  @author Robert Marskar
*/

#ifndef CD_RobinCoefficients_H
#define CD_RobinCoefficients_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract base class for supplying function-based coefficients to the Robin boundary condition classes. Because BaseBCValue is not enough. 
*/
class RobinCoefficients {
public:
  /*!
    @brief Constructor
  */
  RobinCoefficients();

  /*!
    @brief Destructor
  */
  virtual ~RobinCoefficients();

  /*!
    @brief Get a-coefficient
  */
  virtual Real aco(const RealVect a_pos) const = 0;

  /*!
    @brief Get b-coefficient
  */
  virtual Real bco(const RealVect a_pos) const = 0;

  /*!
    @brief Get right-hand side
  */
  virtual Real rhs(const RealVect a_pos) const = 0;
  
protected:
  
};

#include <CD_NamespaceFooter.H>

#endif
