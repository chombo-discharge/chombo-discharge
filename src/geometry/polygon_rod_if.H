/*!
  @file polygon_rod_if.H
  @brief Declaration of a n-sided rod
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POLYGON_ROD_IF_
#define _POLYGON_ROD_IF_

#include <BaseIF.H>

/*!
  @brief Implicit function class describing a hexagonal rod
*/
class polygon_rod_if: public BaseIF {
public:

  /*!
    @brief Base constructor
    @param[in] a_nsides Number of sides of polygon
    @param[in] a_radius Rod radius
    @param[in] a_length Rod length
    @param[in] a_inside Domain on inside (true) or outside (false)
  */
  polygon_rod_if(const int&  a_nsides,
		 const Real& a_radius,
		 const Real& a_length,
		 const Real& a_cornerCurv,
		 const bool& a_inside);
    
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
#endif
