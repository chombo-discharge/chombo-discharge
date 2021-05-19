/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_PerlinPlaneSdf.H
  @brief  Signed distance function for noisy plane
  @author Robert Marskar
*/

#ifndef CD_PerlinPlaneSdf_H
#define CD_PerlinPlaneSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include "CD_NamespaceHeader.H"

/*!
  @brief PerlinPlaneSdf function
*/
class PerlinPlaneSdf : public BaseIF {
public:
  
  /*!
    @brief Constructor
  */
  PerlinPlaneSdf(const RealVect a_normal,
		 const RealVect a_point,
		 const bool     a_inside,
		 const Real     a_noiseAmp,
		 const RealVect a_noiseFreq,
		 const Real     a_persistence,
		 const int      a_octaves,
		 const bool     a_reseed = false);
  /*!
    @brief Copy constructor
  */
  PerlinPlaneSdf(const PerlinPlaneSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~PerlinPlaneSdf();

  /*!
    @brief Level-set function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  /*!
    @brief Plane normal
  */
  RealVect m_normal;

  /*!
    @brief Plane point
  */
  RealVect m_point;

  /*!
    @brief Rod function
  */
  RefCountedPtr<BaseIF> m_plane;

  /*!
    @brief Noise function
  */
  RefCountedPtr<BaseIF> m_perlin;
};

#include <CD_NamespaceFooter.H>

#endif
