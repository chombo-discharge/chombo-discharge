/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxSdf.H
  @brief  Declare an interface for defining an axis-aligned box
  @author Robert Marskar
*/

#ifndef CD_BoxSdf_H
#define CD_BoxSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for defining a two- or three-dimensional box with arbitrary centroid and orientation.
*/
class BoxSdf : public BaseIF
{
public:
  /*!
    @brief Base constructor
    @param[in] a_loCorner Lower left corner
    @param[in] a_hiCorner Upper right corner
    @param[in] a_inside Object on inside or outside
  */
  BoxSdf(const RealVect& a_loCorner, const RealVect& a_hiCorner, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Another BoxSdf
  */
  BoxSdf(const BoxSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~BoxSdf();

  /*!
    @brief Test if object is on inside or outside
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /*!
    @brief IF Factory method
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /*!
    @brief Lower left corner
  */
  RealVect m_loCorner;

  /*!
    @brief Upper right corner
  */
  RealVect m_hiCorner;

  /*!
    @brief Fluid inside box or not. 
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
