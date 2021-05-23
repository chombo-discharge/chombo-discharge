/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_ProfileCylinderIF.H
  @brief   Declaration of a cylinder with surface profile
  @author  Robert Marskar
*/

#ifndef CD_ProfileCylinderIF_H
#define CD_ProfileCylinderIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief ProfileCylinderIF function
*/
class ProfileCylinderIF : public BaseIF {
public:

  ProfileCylinderIF() = delete;

  ProfileCylinderIF(const RealVect  a_endPoint1,
		    const RealVect  a_endPoint2,
		    const Real      a_cylinderRadius,
		    const Real      a_torusMajorRadius,
		    const Real      a_torusMinorRadius,
		    const Real      a_ccDistance,
		    const Real      a_shift,
		    const Real      a_roundingRadius,
		    const int       a_numLeft,
		    const int       a_numRight,
		    const bool      a_fluidInside);
  
  ProfileCylinderIF(const ProfileCylinderIF& a_inputIF);
  
  virtual ~ProfileCylinderIF();
  virtual Real value(const RealVect& a_pos) const override;
  virtual BaseIF* newImplicitFunction() const override;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
