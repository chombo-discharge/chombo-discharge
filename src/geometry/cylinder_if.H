/*!
  @file   cylinder_if.H
  @brief  Declaration of a cylindrical implicit function
  @date   Nov. 2017
  @author Robert Marskar
*/

#ifndef _CYLINDER_IF_
#define _CYLINDER_IF_

#include <BaseIF.H>

namespace ChomboDischarge {

  /*!
    @brief Declaration of a cylinder IF class
  */
  class cylinder_if : public BaseIF {
  public:

    /*!
      @brief Constructor
      @param[in] a_center1 Cylinder center
      @param[in] a_center2 Cylinder center
      @param[in] a_radius  Cylinder radius
      @param[in] a_inside  Domain on inside or outside
    */
    cylinder_if(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_fluidInside);

    /*!
      @brief Copy constructor
      @param[in] a_inputIF Input cylinder
    */
    cylinder_if(const cylinder_if& a_inputIF);

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
      @brief One end poitn
    */
    RealVect m_endPoint1;

    /*!
      @brief The other end point. 
    */
    RealVect m_endPoint2;

    /*!
      @brief Top
    */
    RealVect m_top;

    /*!
      @brief Top
    */
    RealVect m_axis;

    /*!
      @brief Center
    */
    RealVect m_center;

    /*!
      @brief Cylinder length
    */
    Real m_length;

    /*!
      @brief Cylinder radius
    */
    Real m_radius;

    /*!
      @brief Inside or outside domain
    */
    bool m_fluidInside;
  };
}
#endif
