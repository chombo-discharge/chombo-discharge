/*!
  @file   slope_limiters.H
  @brief  Static class with various slope limiters
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _SLOPE_LIMITERS_
#define _SLOPE_LIMITERS_

#include "REAL.H"

/*!
  @brief Class with various slope limiters. Static, so just use slope_limiters::blah()
*/
class slope_limiters {
public:

  /*!
    @brief Koren limiter
  */
  static Real koren(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief minmod limiter
  */
  static Real minmod(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief Superbee limiter
  */
  static Real superbee(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief Van Leer limiter
  */
  static Real van_leer(const Real a_slopeL, const Real a_slopeR);

protected:

  /*!
    @brief Maxmod limiter
  */
  static Real maxmod(const Real a_slopeL, const Real a_slopeR);
};

#endif
