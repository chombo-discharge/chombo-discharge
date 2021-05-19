/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_PerlinSphereSdf.H
  @brief  Declaration of noisy sphere
  @author Robert Marskar
*/

#ifndef CD_PerlinSphereSdf_H
#define CD_PerlinSphereSdf_H

// Chombo includes
#include <BaseIF.H>
#include <SphereIF.H>

// Our includes
#include <CD_PerlinSdf.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief PerlinSphereSdf function
*/
class PerlinSphereSdf : public BaseIF {
public:

  /*!
    @brief No weak construction
  */
  PerlinSphereSdf(){}
  /*!
    @brief Constructor
  */
  PerlinSphereSdf(const Real&     a_rad,
		  const RealVect& a_center,
		  const bool&     a_inside,
		  const Real&     a_noiseAmp,
		  const RealVect& a_noiseFreq,
		  const Real&     a_persistence,
		  const int&      a_octaves,
		  const bool&     a_reseed = false);

  /*!
    @brief Copy constructor
  */
  PerlinSphereSdf(const PerlinSphereSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~PerlinSphereSdf();

  /*!
    @brief Value function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  /*!
    @brief Radius
  */
  Real m_rad;

  /*!
    @brief Center
  */
  RealVect m_center;

  /*!
    @brief Inside/outside
  */
  bool m_inside;
  
  /*!
    @brief Noise function
  */
  RefCountedPtr<BaseIF> m_perlinIF;
};

#include <CD_NamespaceFooter.H>

#endif
