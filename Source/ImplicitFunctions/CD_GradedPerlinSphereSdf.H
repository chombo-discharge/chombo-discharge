/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  @details This class is a sphere with graded noise -- there is only noise on one side of the sphere. 
*/
class GradedPerlinSphereSdf : public PerlinSphereSdf
{
public:
  /*!
    @brief Full constructor
    @param[in] a_rad         Sphere radius
    @param[in] a_center      Sphere center
    @param[in] a_inside      Inside or outside
    @param[in] a_noiseAmp    Noise amplitude
    @param[in] a_noiseFreq   Noise frequency
    @param[in] a_persistence Factor for noise octaves. 
    @param[in] a_octaves     Number of octaves
    @param[in] a_reseed      If true, reseed the RNG. Otherwise, the the original Perlin hash table. 
    @note The persistence argument indicates how we add octaves. With 0.5 we use half the amplitude and twice the frequency for the next octave. 
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
  virtual Real
  value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF*
  newImplicitFunction() const;
};

#include <CD_NamespaceFooter.H>

#endif
