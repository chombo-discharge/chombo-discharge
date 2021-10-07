/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file  CD_HollowCylinderIF.H
  @brief Declaration of a hollow cylinder implicit function.  
  @author Robert Marskar
*/

#ifndef CD_HollowCylinderIF_H
#define CD_HollowCylinderIF_H

// Chombo includes
#include <RealVect.H>
#include <BaseIF.H>
#include <TransformIF.H>
#include <LatheIF.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief This class is used to generate a hollow cylinder with rounded edges. 
*/
class HollowCylinderIF: public BaseIF {
public:

  /*!
    @brief Disallowed weak construction. 
  */
  HollowCylinderIF() = delete;

  /*!
    @brief Full constructor
    @param[in] a_center1        One cylinder end
    @param[in] a_center2        Other cylinder end
    @param[in] a_majorRadius    Outer cylinder radius
    @param[in] a_minorRadius    Inner radius
    @param[in] a_outerCurvature Outer curvature
    @param[in] a_innerCurvature Inner curvature
    @param[in] a_fluidInside Fluid inside or not
  */
  HollowCylinderIF(const RealVect a_center1,
		   const RealVect a_center2,
		   const Real     a_majorRadius,
		   const Real     a_minorRadius,
		   const Real     a_outerCurvature,
		   const Real     a_innerCurvature,		   
		   const bool     a_fluidInside);
  
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder
  */
  HollowCylinderIF(const HollowCylinderIF& a_inputIF);
  
  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const override;

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const override;

protected:

  /*!
    @brief 3D object rotated and translated to its specified position
  */
  RefCountedPtr<BaseIF> m_baseIF;

};

#include <CD_NamespaceFooter.H>

#endif
