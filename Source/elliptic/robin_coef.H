/*!
  @file   robin_coef.H
  @brief  Declaration of a more general class for parametrizing boundary conditions. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _ROBIN_COEF_
#define _ROBIN_COEF_

#include "RealVect.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Abstract base class for supplying function-based coefficients to the Robin boundary condition classes. Because BaseBCValue is not enough. 
*/
class robin_coef {
public:
  /*!
    @brief Constructor
  */
  robin_coef();

  /*!
    @brief Destructor
  */
  virtual ~robin_coef();

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
#include "CD_NamespaceFooter.H"

#endif
