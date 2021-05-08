/*!
  @file   slope_limiters.H
  @brief   class with various slope limiters
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _SLOPE_LIMITERS_
#define _SLOPE_LIMITERS_

#include "REAL.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Class with various slope limiters. , so just use slope_limiters::blah()
*/
namespace slope_limiters {

  /*!
    @brief Koren limiter
  */
  Real koren(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief minmod limiter
  */
  Real minmod(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief Superbee limiter
  */
  Real superbee(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief Van Leer limiter
  */
  Real van_leer(const Real a_slopeL, const Real a_slopeR);

  /*!
    @brief Maxmod limiter
  */
  Real maxmod(const Real a_slopeL, const Real a_slopeR);
};
#include "CD_NamespaceFooter.H"

#endif
