/*!
  @file hollow_cylinder_cross_sec_if.H
  @brief Declaration of an implicit function describing a rectangle with rounded corners
  @date Nov. 2017
  @author Robert Marskar
*/

#ifndef _HOLLOW_CYLINDER_CROSS_SEC_IF_
#define _HOLLOW_CYLINDER_CROSS_SEC_IF_

#include <RealVect.H>
#include <BaseIF.H>

/*!
  @brief The class describes a 2D a rectangle with curved corners (lower left and right). This class is used by HollowCylinderIF to produce a hollow cylinder with rounded bottom edges. 
*/

class hollow_cylinder_cross_sec_if : public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_hiCorner            Upper right corner of the object
    @param[in] a_loCorner            Lower left corner of the object
    @param[in] a_radiusOfCurvature   Radius of curvature of the rounded edges on the object base 
    @param[in] a_inside              Domain on inside or outside
  */
  hollow_cylinder_cross_sec_if(const RealVect& a_hiCorner,
			   const RealVect& a_loCorner,
			   const Real&     a_radiusOfCurvature,
			   const bool&     a_inside);

  

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder cross section
  */
  hollow_cylinder_cross_sec_if(const hollow_cylinder_cross_sec_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~hollow_cylinder_cross_sec_if();

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
  /*!
    @brief Upper right corner of the cross section 
  */
  RealVect m_hiCorner;

  /*!
    @brief Lower left corner of the cross section 
  */
  RealVect m_loCorner;
  
  /*!
    @brief Radius of curvature of the rounded edges on the cylinder base 
  */ 
  Real m_radOfCurv;
  
  /*!
    @brief Inside or outside domain
  */   
  bool m_inside;

};

#endif
