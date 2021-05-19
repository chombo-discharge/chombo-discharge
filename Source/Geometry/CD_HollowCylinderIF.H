/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file  CD_HollowCylinderIF.H
  @brief Declaration of a hollow cylinder implicit function.  
  @author Robert Marskar
*/

#ifndef CD_HollowCylinderIF_H
#define CD_HollowCylinderIF_H

// Chombo includes
#include <RealVect.H>
#include <BaseIF.H>
#include <TransformIF.H>
#include <LatheIF.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief This class is used to generate a hollow cylinder. The edges on the base of the cylinder are rounded, while the edges on the cylinder top are sharp. 
*/
class HollowCylinderIF: public BaseIF {
public:

  /*!
    @brief Oh no you dont!
  */
  HollowCylinderIF() = delete;

  /*!
    @brief ctor
  */
  HollowCylinderIF(const RealVect a_center1,
		   const RealVect a_center2,
		   const Real     a_majorRadius,
		   const Real     a_minorRadius,
		   const Real     a_curv,
		   const bool     a_fluidInside);
  
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input hollow cylinder
  */
  HollowCylinderIF(const HollowCylinderIF& a_inputIF);
  
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

#include <CD_NamespaceFooter.H>

#endif
