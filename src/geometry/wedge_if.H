/*!
  @file wedge_if.H
  @brief Declaration of a n-sided rod
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _WEDGE_IF_
#define _WEDGE_IF_
#include <BaseIF.H>

namespace ChomboDischarge {

  /*!
    @brief Implicit function class describing a hexagonal rod
  */
  class wedge_if: public BaseIF {
  public:

    /*!
      @brief Base constructor
      @param[in] a_dir Direction, 0=x,1=y,2=z
      @param[in] a_angle Opening angle
      @param[in] a_curv Corner curvature
      @param[in] a_point Center of circle used for curvature
    */
    wedge_if(const int      a_dir,
	     const Real     a_angle,
	     const Real     a_curv,
	     const RealVect a_point,
	     const bool     a_inside);
    
    /*!
      @brief Copy constructor
      @param[in] a_inputIF Another polygon rod
    */
    wedge_if(const wedge_if& a_inputIF);

    /*!
      @brief Destructor
    */
    virtual ~wedge_if();

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
}
#endif
