/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_GradedPerlinSphereSdf.H
  @brief  Declaration of a graded noisy sphere
  @author Robert Marskar
*/

#ifndef CD_GradedPerlinSphereSdf_H
#define CD_GradedPerlinSphereSdf_H

// Our includes
#include <CD_PerlinSphereSdf.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief GradedPerlinSphereSdf function
*/
class GradedPerlinSphereSdf : public PerlinSphereSdf {
public:  
  /*!
    @brief Constructor
  */
  GradedPerlinSphereSdf(const Real&     a_rad,
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

  GradedPerlinSphereSdf(const GradedPerlinSphereSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~GradedPerlinSphereSdf();

  /*!
    @brief Value function
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF* newImplicitFunction() const;
};

#include <CD_NamespaceFooter.H>

#endif
