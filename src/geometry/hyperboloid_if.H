/*!
  @file hyperboloid_if.H
  @brief Declaration of a single-sheet hyperboloidal surface from a two-sheet hyperboloid
  @date Nov. 2017
  @author Robert Marskar
*/

#ifndef _HYPERBOLOIDIF_
#define _HYPERBOLOIDIF_

#include <RealVect.H> 
#include <BaseIF.H>

/*!
  @brief Declaration of a two-sheet hyperboloid. 
  @details The implicit function used for the single sheet is
  \f[
  \frac{z}{c} - \sqrt{1 + \frac{x^2}{a^2} + \frac{y^2}{b^2}} = 0
  \f]
*/
class hyperboloid_if : public BaseIF{
public:

  /*!
    @brief Constructor
    @param[in] a_radii Hyperboloid radii
    @param[in] a_center Hyperboloid center
    @param[in] a_inside Domain on inside or outside
  */
  hyperboloid_if(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input one-sheeted hyperboloid
  */
  hyperboloid_if(const hyperboloid_if& a_inputIF);

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const; 

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const;

protected:
  /*!
    @brief Hyperboloid radii
  */
  RealVect m_radii;

  /*!
    @brief Hyperboloid radii squared
  */
  RealVect m_radii2;

  /*!
    @brief Holder for sign convention
  */
  RealVect m_sign;

  /*!
    @brief Hyperboloid center
  */
  RealVect m_center;

  /*!
    @brief Inside or outside domain
  */
  bool m_inside;
};
#endif
