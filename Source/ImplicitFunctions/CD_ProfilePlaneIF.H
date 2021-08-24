/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_ProfilePlaneIF.H
  @brief  Declaration of an implicit function which is a plane with "holes"
  @author Robert Marskar
*/

#ifndef CD_ProfilePlaneIF_H
#define CD_ProfilePlaneIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Implicit function which represents an infinite half-plane with "holes". 
  @note This is currently limited to 2D. 
*/
class ProfilePlaneIF : public BaseIF {
public:

  /*!
    @brief Disallowed weak construction
  */
  ProfilePlaneIF() = delete;

  /*!
    @brief Full constructor.
    @param[in] a_point       Point on the plane
    @param[in] a_width       Box width
    @param[in] a_impFunc     Implicit function to extract from the box
    @param[in] a_numLeft     Number of holes to the left (of the centerpoint)
    @param[in] a_numRight    Number of holes to the right (of the centerpoint)
    @param[in] a_ccDist      Distance between holes
    @param[in] a_xShift      Profile shift to the left/right (relative to centerponit)
    @param[in] a_yShift      Profile shift up/down (relative to centerponit)
    @param[in] a_curv        Corner curvatures
    @param[in] a_fluidInside Fluid on inside or not
  */
  ProfilePlaneIF(const RealVect  a_point,
		 const Real      a_width,
		 const BaseIF*   a_impFunc,
		 const int       a_numLeft,
		 const int       a_numRight,
		 const Real      a_ccDist,
		 const Real      a_xShift,
		 const Real      a_yShift,
		 const Real      a_curv,
		 const bool      a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Other plane
  */
  ProfilePlaneIF(const ProfilePlaneIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~ProfilePlaneIF();

  /*!
    @brief Level-set function
    @param[in] a_pos Position
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;

  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
