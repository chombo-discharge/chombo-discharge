/*!
  @file hyperboloid_one_if.H
  @brief Declaration of a one-sheet hyperboloidal shape
  @date Nov. 2017
  @author Robert Marskar
*/

#ifndef _HYPERBOLOID_ONE_IF_
#define _HYPERBOLOID_ONE_IF_

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Declaration of a one-sheet hyperboloid. 
  @details The implicit function used for the sheet is 
  \f[
  \frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac{z^2}{c^2} = 1.
  \f]
*/
class hyperboloid_one_if : public BaseIF{
public:

  /*!
    @brief Constructor
    @param[in] a_radii Hyperboloid radii
    @param[in] a_center Hyperboloid center
    @param[in] a_inside Domain on inside or outside
  */
  hyperboloid_one_if(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input one-sheeted hyperboloid
  */
  hyperboloid_one_if(const hyperboloid_one_if& a_inputIF);

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
#include "CD_NamespaceFooter.H"
#endif
