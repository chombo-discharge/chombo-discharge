/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  @brief Noisy sphere geometry (with Perlin noise).
*/
class PerlinSphereSdf : public BaseIF
{
public:
  /*!
    @brief No weak construction
  */
  PerlinSphereSdf()
  {}

  /*!
    @brief Constructor
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
    @param[in] a_inputIF Other perlin sphere
  */
  PerlinSphereSdf(const PerlinSphereSdf& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~PerlinSphereSdf();

  /*!
    @brief Value function. Returns distance to object
    @param[in] a_pos Position
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /*!
    @brief Factory function
  */
  virtual BaseIF*
  newImplicitFunction() const;

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
