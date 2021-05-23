/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PolygonRodIF.H
  @brief  Declaration of a n-sided rod
  @author Robert Marskar
*/

#ifndef CD_PolygonRodIF_H
#define CD_PolygonRodIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Implicit function class describing a hexagonal rod
*/
class PolygonRodIF: public BaseIF {
public:

  /*!
    @brief Base constructor
    @param[in] a_endPoint1  One endpoint
    @param[in] a_endPoint2  Other endpoint
    @param[in] a_length     Rod length
    @param[in] a_cornerCurv Rounding radius on corners
    @param[in] a_numSides   Number of sides
    @param[in] a_inside Domain on inside (true) or outside (false)
  */
  PolygonRodIF(const RealVect a_center1,
	       const RealVect a_center2,
	       const Real     a_radius,
	       const Real     a_cornerCurv,
	       const int      a_numSides,
	       const bool     a_fluidInside);

    
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Another polygon rod
  */
  PolygonRodIF(const PolygonRodIF& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~PolygonRodIF();

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
