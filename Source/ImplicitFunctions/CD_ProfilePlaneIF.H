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
  @brief ProfilePlaneIF function
*/
class ProfilePlaneIF : public BaseIF {
public:

  /*!
    @brief Ctor. 
  */
  ProfilePlaneIF(const RealVect  a_point,
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
  ProfilePlaneIF(const ProfilePlaneIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~ProfilePlaneIF();

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

#include <CD_NamespaceFooter.H>

#endif
