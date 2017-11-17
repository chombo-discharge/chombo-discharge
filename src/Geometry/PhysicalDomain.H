/*!
  @file PhysicalDomain.H
  @brief Physical domain class. Very lightweight stuff.
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PHYSICALDOMAIN_
#define _PHYSICALDOMAIN_

#include <RealVect.H>

class PhysicalDomain {
public:
  /*!
    @brief Constructor
  */
  PhysicalDomain();

  /*!
    @brief Full constructor
  */
  PhysicalDomain(const RealVect& a_probLo, const RealVect& a_probHi);

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
