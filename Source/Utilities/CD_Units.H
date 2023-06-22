/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Units.H
  @brief  Declaration of various useful units
  @author Robert Marskar
*/

#ifndef CD_Units_H
#define CD_Units_H

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Various SI units
*/
namespace Units {

  /*!
    @brief Permittivity of free space
  */
  constexpr Real eps0 = 8.854187817E-12;

  /*!
    @brief Elementary charge 
  */
  constexpr Real Qe = 1.6021766208E-19;

  /*!
    @brief Speed of light
  */
  constexpr Real c = 299792458;

  /*!
    @brief Boltzmann constant
  */
  constexpr Real kb = 1.38084852E-23;

  /*!
    @brief Towsend unit
  */
  constexpr Real Td = 1.E-21;

  /*!
    @brief Electron volt
  */
  constexpr Real eV = 1.6021766208E-19;

  /*!
    @brief Avogadro number
  */
  constexpr Real Na = 6.022140857E23;

  /*!
    @brief Universal gas constant
  */
  constexpr Real R = 8.3144598;

  /*!
    @brief Conversion from atm to pascal
  */
  constexpr Real atm2pascal = 101325;

  /*!
    @brief pi
  */
  constexpr Real pi = 3.14159265359;

  /*!
    @brief Electron mass
  */
  constexpr Real me = 9.109383E-31;
} // namespace Units

#include <CD_NamespaceFooter.H>

#endif
