/*!
  @file   rounded_cylinder_if.H
  @brief  Declaration of a cylindrical implicit function with rounded edges
  @date   Feb. 2021
  @author Robert Marskar
*/

#ifndef _ROUNDED_CYLINDER_IF_
#define _ROUNDED_CYLINDER_IF_

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Declaration of a cylinder IF class
*/
class rounded_cylinder_if : public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_center1 Cylinder center
    @param[in] a_center2 Cylinder center
    @param[in] a_radius  Cylinder radius
    @param[in] a_curv    Radius of curvature
    @param[in] a_inside  Domain on inside or outside
  */
  rounded_cylinder_if(const RealVect a_center1, const RealVect a_center2, const Real a_radius, const Real a_curv, const bool a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  rounded_cylinder_if(const rounded_cylinder_if& a_inputIF);

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const override; 

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const override;

protected:

  /*!
    @brief Cylinder center
  */
  RealVect m_center1;

  /*!
    @brief Cylinder center
  */
  RealVect m_center2;

  /*!
    @brief Cylinder length
  */
  Real m_length;

  /*!
    @brief Cylinder radius
  */
  Real m_radius;

  /*!
    @brief Radius of curvature
  */
  Real m_curv;

  /*!
    @brief Inside or outside domain
  */
  bool m_fluidInside;

  /*!
    @brief Implicit function
  */
  RefCountedPtr<BaseIF> m_baseif;

  /*!
    @brief Make the implicit function
  */
  void makeBaseIF();

  /*!
    @brief Make BaseIF
  */
#if CH_SPACEDIM==2
  BaseIF* makeBaseIF2D();
#elif CH_SPACEDIM==3
  BaseIF* makeBaseIF3D();
#endif
};
#endif
#include "CD_NamespaceFooter.H"
