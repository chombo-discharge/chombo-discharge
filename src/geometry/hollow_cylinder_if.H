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

#include "CD_NamespaceHeader.H"
  
/*!
  @brief This class is used to generate a hollow cylinder. The edges on the base of the cylinder are rounded, while the edges on the cylinder top are sharp. 
*/
class hollow_cylinder_if: public BaseIF {
public:

  /*!
    @brief Oh no you dont!
  */
  hollow_cylinder_if() = delete;

  /*!
    @brief ctor
  */
  hollow_cylinder_if(const RealVect a_center1,
		     const RealVect a_center2,
		     const Real     a_majorRadius,
		     const Real     a_minorRadius,
		     const Real     a_curv,
		     const bool     a_fluidInside);
  
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder
  */
  hollow_cylinder_if(const hollow_cylinder_if& a_inputIF);
  
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
    @brief 3D object rotated and translated to its specified position
  */
  RefCountedPtr<BaseIF> m_baseif;

};
#include "CD_NamespaceFooter.H"

#endif
