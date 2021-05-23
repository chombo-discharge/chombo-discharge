/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_WedgeIF.H
  @brief  Wedge implicit function
  @author Robert Marskar
*/

#ifndef CD_WedgeIF_H
#define CD_WedgeIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Implicit function class describing a hexagonal rod
*/
class WedgeIF: public BaseIF {
public:

  /*!
    @brief Base constructor
    @param[in] a_dir Direction, 0=x,1=y,2=z
    @param[in] a_angle Opening angle
    @param[in] a_curv Corner curvature
    @param[in] a_point Center of circle used for curvature
  */
  WedgeIF(const int      a_dir,
	   const Real     a_angle,
	   const Real     a_curv,
	   const RealVect a_point,
	   const bool     a_inside);
    
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Another polygon rod
  */
  WedgeIF(const WedgeIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~WedgeIF();

  /*!
    @brief Test if object is on inside or outside
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief IF Factory method
  */
  virtual BaseIF* newImplicitFunction() const;

protected:

  /*!
    @brief Base function
  */
  RefCountedPtr<BaseIF> m_baseif;
  
};

#include <CD_NamespaceFooter.H>

#endif
