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

#include "CD_NamespaceHeader.H"

/*!
  @brief profile_plane_if function
*/
class profile_plane_if : public BaseIF {
public:

  /*!
    @brief Ctor. 
  */
  profile_plane_if(const RealVect  a_point,
		   const Real      a_width,
		   const BaseIF*   a_impFunc,
		   const int       a_num_left,
		   const int       a_num_right,
		   const Real      a_ccDist,
		   const Real      a_xShift,
		   const Real      a_yShift,
		   const Real      a_curv,
		   const bool      a_fluidInside);

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

  bool m_fluidInside;
};
#include "CD_NamespaceFooter.H"
#endif
