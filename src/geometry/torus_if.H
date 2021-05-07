/*!
  @file   torus_if.H
  @brief  New TorusIF function, which is now also signed distance
  @date   Sep. 2017
  @author Robert Marskar
*/

#ifndef _NEW_TORUS_IF_
#define _NEW_TORUS_IF_

#include <BaseIF.H>

namespace ChomboDischarge {

  /*!
    @brief Declaration of a torus signed distance function. The orientation is along z
  */
  class torus_if : public BaseIF {
  public:

    torus_if() = delete;

    /*!
      @brief Constructor
      @param[in] a_center      Center
      @param[in] a_majorRadius Major radius
      @param[in] a_minorRadius Minor radius
      @param[in] a_fluidInside Fluid on inside of or not
    */
    torus_if(const RealVect a_center, const Real a_majorRadius, const Real a_minorRadius, const bool a_fluidInside);

    /*!
      @brief Copy constructor
      @param[in] a_inputIF Input cylinder
    */
  
    torus_if(const torus_if& a_inputIF);

    /*!
      @brief Destructor
    */
    virtual ~torus_if();

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
      @brief Sphere center
    */
    RealVect m_center;

    /*!
      @brief Cylinder radius
    */
    Real m_majorRadius;

    /*!
      @brief Cylinder radius
    */
    Real m_minorRadius;
  
    /*!
      @brief Inside or outside domain
    */
    bool m_fluidInside;
  };
}
#endif
