/*!
  @file   real_box.H
  @brief  Declaration of a rectangular box
  @author Robert Marskar
  @date   June 2018
*/

#ifndef _REAL_BOX_
#define _REAL_BOX_

#include "RealVect.H"

/*!
  @brief Rectangular box
*/
class real_box{
public:

  /*!
    @brief Constructor
  */
  real_box(){
    m_lo = RealVect::Zero;
    m_hi = RealVect::Zero;
  }

  /*!
    @brief Constructor
  */
  real_box(const RealVect a_lo, const RealVect a_hi){
    m_lo = a_lo;
    m_hi = a_hi;
  }

  /*!
    @brief Destructor
  */
  ~real_box(){

  }

  /*!
    @brief Get lo corner
  */
  RealVect get_lo() const {
    return m_lo;
  }

  /*!
    @brief Get hi corner
  */
  RealVect get_hi() const {
    return m_hi;
}
  
protected:

  RealVect m_lo;
  RealVect m_hi;

};

#endif
