/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_PerlinRodSdf.H
  @brief  Declaration of a rod with a rounded, noisy tip
  @author Robert Marskar
*/

#ifndef CD_PerlinRodSdf_H
#define CD_PerlinRodSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_CylinderSdf.H>
#include <CD_GradedPerlinSphereSdf.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief PerlinRodSdf function
*/
class PerlinRodSdf : public PerlinSphereSdf {
public:  
  /*!
    @brief Constructor
  */
  PerlinRodSdf(const Real&     a_rad,
		const RealVect& a_center1,
		const RealVect& a_center2,
		const bool&     a_inside,
		const Real&     a_noiseAmp,
		const RealVect& a_noiseFreq,
		const Real&     a_persistence,
		const int&      a_octaves,
		const bool&     a_reseed = false);

  /*!
    @brief Copy constructor
  */

  PerlinRodSdf(const PerlinRodSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~PerlinRodSdf();

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
    @brief Basis rod, used for copying
  */
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
