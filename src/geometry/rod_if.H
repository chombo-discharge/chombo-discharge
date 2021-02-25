/*!
  @file rod_if.H
  @brief Declaration of a cylinder with spheres at both ends
  @date Nov. 2017
  @author Robert marskar
*/

#ifndef _ROD_IF_
#define _ROD_IF_

#include <BaseIF.H>
#include <IntersectionIF.H>

#include <ParmParse.H>

/*!
  @brief Declaration of a cylinder with rounded caps
*/
class rod_if : public BaseIF{
public:

  /*!
    @brief Constructor
    @param[in] a_center1 Cylinder center (at center of sphere)
    @param[in] a_center2 Cylinder center (at center of sphere)
    @param[in] a_radius  Cylinder and sphere radius
    @param[in] a_fluidInside Domain on inside or outside
  */

  rod_if(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  
  rod_if(const rod_if& a_inputIF);

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

  RefCountedPtr<BaseIF> m_baseif;
};
#endif
