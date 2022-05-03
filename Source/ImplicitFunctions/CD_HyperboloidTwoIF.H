/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_HyperboloidTwoIF.H
  @brief  Declaration of a two-sheet hyperboloidal implicit function object
  @author Robert Marskar
*/

#ifndef CD_HyperboloidTwoIF_H
#define CD_HyperboloidTwoIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a two-sheet hyperboloid. 
  @details The implicit function used for the sheets is 
  \f[
  \frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac{z^2}{c^2} = -1
  \f]
*/
class HyperboloidTwoIF : public BaseIF
{
public:
  /*!
    @brief Constructor
    @param[in] a_radii Hyperboloid radii
    @param[in] a_center Hyperboloid center
    @param[in] a_inside Domain on inside or outside
  */
  HyperboloidTwoIF(const RealVect& a_radii, const RealVect& a_center, const bool& a_inside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input one-sheeted hyperboloid
  */
  HyperboloidTwoIF(const HyperboloidTwoIF& a_inputIF);

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real
  value(const RealVect& a_point) const;

  /*!
    @brief IF factory method
  */
  virtual BaseIF*
  newImplicitFunction() const;

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

#include <CD_NamespaceFooter.H>

#endif
