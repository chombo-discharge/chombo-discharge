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
class RoundedBoxIF : public BaseIF {
public:

  RoundedBoxIF() = delete;
  
  /*!
    @brief Constructor
    @param[in] a_loCorner            Upper right corner of the object
    @param[in] a_hiCorner            Lower left corner of the object
    @param[in] a_curvature           Radius of curvature on edges and corners
    @param[in] a_fluidInside         Fluid on inside (true) or outside (false)
  */
  RoundedBoxIF(const RealVect a_loCorner,
	       const RealVect a_hiCorner,
	       const Real     a_curvature,
	       const bool     a_fluidInside);


  /*!
    @brief Copy constructor
  */
  RoundedBoxIF(const RoundedBoxIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~RoundedBoxIF();

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

  bool m_fluidInside;
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
