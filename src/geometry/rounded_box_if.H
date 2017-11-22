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

  /*!
    @brief Constructor
    @param[in] a_loCorner            Upper right corner of the object
    @param[in] a_hiCorner            Lower left corner of the object
    @param[in] a_curvature           Radius of curvature on edges and corners
    @param[in] a_inside              Domain on inside or outside
  */
  rounded_box_if(const RealVect& a_loCorner,
		 const RealVect& a_hiCorner,
		 const Real&     a_curvature,
		 const bool&     a_inside);


  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder cross section
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

#if CH_SPACEDIM == 2
  /*!
    @brief Build 2D box
  */
  virtual void buildBox2D();
#elif CH_SPACEDIM == 3
  /*!
    @brief Build 3D box
  */
  virtual void buildBox3D();
#endif
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
  Real m_curv;
  
  /*!
    @brief Inside or outside domain
  */   
  bool m_inside;

  /*!
    @brief Actual rounded box, built up by a bunch of boolean operations
  */
  RefCountedPtr<BaseIF> m_baseif;
};

#endif
