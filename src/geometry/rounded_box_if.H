/*!
  @file rounded_box_if.H
  @brief Declaration of an implicit function describing a rectangle with rounded corners (and edges in 3D)
  @date Nov. 2017
  @author Robert Marskar
*/

#ifndef _ROUNDED_BOX_IF_
#define _ROUNDED_BOX_IF_

#include <RealVect.H>
#include <BaseIF.H>

/*!
  @brief A box with rounded corners and edges
*/
class rounded_box_if : public BaseIF {
public:

  rounded_box_if() = delete;
  
  /*!
    @brief Constructor
    @param[in] a_loCorner            Upper right corner of the object
    @param[in] a_hiCorner            Lower left corner of the object
    @param[in] a_curvature           Radius of curvature on edges and corners
    @param[in] a_fluidInside         Fluid on inside (true) or outside (false)
  */
  rounded_box_if(const RealVect a_loCorner,
		 const RealVect a_hiCorner,
		 const Real     a_curvature,
		 const bool     a_fluidInside);


  /*!
    @brief Copy constructor
  */
  rounded_box_if(const rounded_box_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~rounded_box_if();

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const;

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  bool m_fluidInside;
  RefCountedPtr<BaseIF> m_baseif;
};

#endif
