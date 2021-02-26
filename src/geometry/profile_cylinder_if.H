/*!
  @file    profile_cylinder_if.H
  @brief   Declaration of a cylinder with surface profile
  @author  Robert Marskar
*/

#ifndef _PROFILE_CYLINDER_IF_
#define _PROFILE_CYLINDER_IF_

#include <BaseIF.H>

/*!
  @brief profile_cylinder_if function
*/
class profile_cylinder_if : public BaseIF {
public:

  profile_cylinder_if() = delete;

  profile_cylinder_if(const RealVect  a_endPoint1,
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
  
  profile_cylinder_if(const profile_cylinder_if& a_inputIF);
  
  virtual ~profile_cylinder_if();
  virtual Real value(const RealVect& a_pos) const override;
  virtual BaseIF* newImplicitFunction() const override;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;
};
#endif
