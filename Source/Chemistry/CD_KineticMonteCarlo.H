/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KineticMonteCarlo.H
  @brief  Class for running Kinetic Monte Carlo functionality. 
  @author Robert Marskar
*/

#ifndef CD_KineticMonteCarlo_H
#define CD_KineticMonteCarlo_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for running Kinetic Monte-Carlo simulations. 
  @details The template parameters indicate the underlying reaction type
*/
template <typename R>
class KineticMonteCarlo
{
public:
  /*!
    @brief Default constructor
  */
  KineticMonteCarlo() = default;

  /*!
    @brief Destructor
  */
  virtual ~KineticMonteCarlo();

protected:
};

#include <CD_NamespaceFooter.H>

#endif
