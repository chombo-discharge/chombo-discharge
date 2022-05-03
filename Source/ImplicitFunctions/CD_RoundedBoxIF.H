/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoundedBoxIF.H
  @brief  Declaration of an implicit function describing a rectangle with rounded corners (and edges in 3D)
  @author Robert Marskar
*/

#ifndef CD_RoundedBoxIF_H
#define CD_RoundedBoxIF_H

// Chombo includes
#include <RealVect.H>
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief A box with rounded corners and edges
*/
class RoundedBoxIF : public BaseIF
{
public:
  /*!
    @brief Disallowed weak construction
  */
  RoundedBoxIF() = delete;

  /*!
    @brief Constructor
    @param[in] a_loCorner    Upper right corner of the object
    @param[in] a_hiCorner    Lower left corner of the object
    @param[in] a_curvature   Radius of curvature on edges and corners
    @param[in] a_fluidInside Fluid on inside (true) or outside (false)
  */
  RoundedBoxIF(const RealVect a_loCorner, const RealVect a_hiCorner, const Real a_curvature, const bool a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Other box
  */
  RoundedBoxIF(const RoundedBoxIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~RoundedBoxIF();

  /*!
    @brief Get distance to object. 
    @param[in] a_point Physical position. 
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /*!
    @brief IF factory method
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief Fluid on inside or not
  */
  bool m_fluidInside;

  /*!
    @brief Implicit function (which is union between planes)
  */
  RefCountedPtr<BaseIF> m_baseIF;
};

#include <CD_NamespaceFooter.H>

#endif
