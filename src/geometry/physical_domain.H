/*!
  @file physical_domain.H
  @brief Physical domain class. Very lightweight stuff.
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PHYSICAL_DOMAIN_
#define _PHYSICAL_DOMAIN_

#include <RealVect.H>

class physical_domain {
public:
  /*!
    @brief Constructor
  */
  physical_domain();

  /*!
    @brief Full constructor
  */
  physical_domain(const RealVect& a_probLo, const RealVect& a_probHi);

  /*!
    @brief Define function
   */
  void define(const RealVect& a_probLo, const RealVect& a_probHi);

  /*!
    @brief Get origin
  */
  const RealVect& get_prob_lo() const;

  /*!
    @brief Get high end
  */
  const RealVect& get_prob_hi() const;
  
protected:

  /*!
    @brief Low end
  */
  RealVect m_probLo;

  /*!
    @brief Hi end
  */
  RealVect m_probHi;
};
#endif
