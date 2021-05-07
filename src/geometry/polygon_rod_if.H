/*!
  @file polygon_rod_if.H
  @brief Declaration of a n-sided rod
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POLYGON_ROD_IF_
#define _POLYGON_ROD_IF_


#include <BaseIF.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Implicit function class describing a hexagonal rod
*/
class polygon_rod_if: public BaseIF {
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
  polygon_rod_if(const RealVect a_center1,
		 const RealVect a_center2,
		 const Real     a_radius,
		 const Real     a_cornerCurv,
		 const int      a_numSides,
		 const bool     a_fluidInside);

    
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Another polygon rod
  */
  polygon_rod_if(const polygon_rod_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~polygon_rod_if();

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
#include "CD_NamespaceFooter.H"
#endif
