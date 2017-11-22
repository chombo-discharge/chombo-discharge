/*!
  @file hollow_cylinder_if.H
  @brief Declaration of a hollow cylinder implicit function.  
  @date Nov. 2017
  @author Robert Marskar
*/

#ifndef _HOLLOW_CYLINDER_IF_
#define _HOLLOW_CYLINDER_IF_

#include <RealVect.H>
#include <BaseIF.H>
#include <TransformIF.H>
#include <LatheIF.H>

/*!
  @brief This class is used to generate a hollow cylinder. The edges on the base of the cylinder are rounded, while the edges on the cylinder top are sharp. 
*/
class hollow_cylinder_if: public BaseIF {
public:
  
  /*!
    @brief Constructor
    @param[in] a_innerRadius         Inner cylinder radius
    @param[in] a_outerRadius         Outer Cylinder radius
    @param[in] a_height              Cylinder height
    @param[in] a_radiusOfCurvature   Radius of curvature of the rounded edges on the cylinder base 
    @param[in] a_center              Coordinates of the center of the cylinder
    @param[in] a_inside              Domain on inside or outside
    @param[in] a_dir                 Default orientation
  */
  hollow_cylinder_if(const Real&     a_innerRadius,
		     const Real&     a_outerRadius,
		     const Real&     a_height,
		     const Real&     a_radiusOfCurvature,
		     const RealVect& a_center,
		     const bool&     a_inside,
		     const int&      a_dir = 2);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder
  */
  hollow_cylinder_if(const hollow_cylinder_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~hollow_cylinder_if();

  /*!
    @brief Uses the class HollowCylCrossSec to produce the 3D object and stores it as m_object.
  */
  void buildHollowCylinder();
  
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
    @brief Inner cylinder radius 
  */
  Real m_innerRadius;

  /*!
    @brief Outer cylinder radius 
  */
  Real m_outerRadius;

  /*!
    @brief Cylinder height
  */  
  Real m_height;
  
  /*!
    @brief Radius of curvature of the rounded edges on the cylinder base 
  */ 
  Real m_radiusOfCurvature;
  
  /*!
    @brief Coordinates of the center of the cylinder 
  */   
  RealVect m_center;

  /*!
    @brief Inside or outside domain
  */   
  bool m_inside;

  /*!
    @brief 3D object rotated and translated to its specified position
  */
  TransformIF* m_object;

  /*!
    @brief Orientation of object
  */
  int m_dir;
  
private:

  hollow_cylinder_if(){
    MayDay::Abort("hollow_cylinder_if uses strong construction");
  }

  void operator=(const hollow_cylinder_if& a_inputIF){
    MayDay::Abort("hollow_cylinder_if doesn't allow assignment");
  }
};

#endif
