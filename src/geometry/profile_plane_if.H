/*!
  @file profile_plane_if.H
  @brief Declaration of a noise plane level-set function
  @details This class is essentially a PlaneIF with profile_if noise
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PROFILE_PLANE_IF_
#define _PROFILE_PLANE_IF_

#include <BaseIF.H>

/*!
  @brief profile_plane_if function
*/
class profile_plane_if : public BaseIF {
public:  
  /*!
    @brief Constructor
  */
  profile_plane_if(const RealVect  a_point,
		   const Real      a_width,
		   const int       a_num_left,
		   const int       a_num_right,
		   const Real      a_rad,
		   const Real      a_offset,
		   const Real      a_shift,
		   const Real      a_dist,
		   const Real      a_curv,
		   const bool      a_inside);
  /*!
    @brief Copy constructor
  */
  profile_plane_if(const profile_plane_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~profile_plane_if();

  /*!
    @brief Level-set function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  RefCountedPtr<BaseIF> m_baseif;
};
#endif
