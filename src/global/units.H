/*!
  @file units.H
  @brief Declaration of static natural constants in SI units
  @author Robert Marskar
  @date Jan 2018
*/

#ifndef _UNITS_
#define _UNITS_

namespace ChomboDischarge {
  /*!
    @brief Various SI units
  */
  namespace units {

    /*!
      @brief Permittivity of free space
    */
    const static Real s_eps0 = 8.854187817E-12;

    /*!
      @brief Elementary charge 
    */
    const static Real s_Qe   = 1.6021766208E-19;

    /*!
      @brief Speed of light
    */
    const static Real s_c0   = 299792458;

    /*!
      @brief Boltzmann constant
    */
    const static Real s_kb   = 1.38084852E-23;

    /*!
      @brief Towsend unit
    */
    const static Real s_Td   = 1.E-21;

    /*!
      @brief Electron volt
    */
    const static Real s_eV   = s_Qe;

    /*!
      @brief Avogadro number
    */
    const static Real s_Na   = 6.022140857E23;

    /*!
      @brief Universal gas constant
    */
    const static Real s_R    = 8.3144598;

    /*!
      @brief Conversion from atm to pascal
    */
    const static Real s_atm2pascal  = 101325;

    /*!
      @brief pi
    */
    const static Real s_pi = 3.14159265359;

    /*!
      @brief Electron mass
    */
    const static Real s_me = 9.109383E-31;
  }
}
#endif
