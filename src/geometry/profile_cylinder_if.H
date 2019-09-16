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
  
  profile_cylinder_if(const RealVect  a_center,
		      const Real      a_length,
		      const Real      a_radius,
		      const int       a_num_left,
		      const int       a_num_right,
		      const Real      a_rad,
		      const Real      a_offset,
		      const Real      a_shift,
		      const Real      a_dist,
		      const Real      a_curv,
		      const bool      a_inside);
  profile_cylinder_if(const profile_cylinder_if& a_inputIF);
  
  virtual ~profile_cylinder_if();
  virtual Real value(const RealVect& a_pos) const;
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;
};
#endif
