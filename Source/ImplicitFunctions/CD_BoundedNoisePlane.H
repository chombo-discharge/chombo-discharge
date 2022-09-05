/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoundedNoisePlane.H
  @brief  Declaration of a signed distance function for a noisy plane
  @author Robert Marskar
*/

#ifndef CD_BoundedNoisePlane_H
#define CD_BoundedNoisePlane_H

// Chombo includes
#include <PlaneIF.H>
#include <BaseIF.H>

// Our includes
#include <CD_PerlinSdf.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Infinite plane with spatially bounded surface noise (Perlin noise). 
  @details The surface noise is bounded by a smoothed boxcar function between the input arguments. The boxcar function is given by
  b(x) = h(x-a) - h(x-b) where h(x) = 1/(1 + exp(-2*k*x)) where k determines the falloff.
*/
class BoundedNoisePlane : public BaseIF
{
public:
  /*!
    @brief Full constructor
    @param[in] a_normal      Plane normal vector
    @param[in] a_point       Point on the plane
    @param[in] a_clampLo     Low clamped corner. 
    @param[in] a_clampHi     High clamped corner.
    @param[in] a_clampK      K-value for clamping.
    @param[in] a_noiseAmp    Noise amplitude
    @param[in] a_noiseFreq   Noise frequency
    @param[in] a_persistence Factor for noise octaves. 
    @param[in] a_octaves     Number of octaves
    @param[in] a_reseed      If true, reseed the RNG. Otherwise, the the original Perlin hash table. 
    @note The persistence argument indicates how we add octaves. With 0.5 we use half the amplitude and twice the frequency for the next octave. 
  */
  BoundedNoisePlane(const int      a_normal,
                    const RealVect a_point,
                    const RealVect a_clampLo,
                    const RealVect a_clampHi,
                    const Real     a_clampK,
                    const Real     a_noiseAmp,
                    const RealVect a_noiseFreq,
                    const Real     a_persistence,
                    const int      a_octaves,
                    const bool     a_reseed = false);
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Other plane
  */
  BoundedNoisePlane(const BoundedNoisePlane& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~BoundedNoisePlane();

  /*!
    @brief Value function
    @param[in] a_pos Position
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
    @return Returns new plane with same parameters. 
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /*!
    @brief Plane direction
  */
  int m_normal;

  /*!
    @brief Plane point
  */
  RealVect m_point;

  /*!
    @brief Low clamped value
  */
  RealVect m_clampLo;

  /*!
    @brief Low clamped value
  */
  RealVect m_clampHi;

  /*!
    @brief Transition length
  */
  Real m_clampK;

  /*!
    @brief Maximum noise amplitude
  */
  Real m_maxAmp;

  /*!
    @brief BaseIF function for a plane without noise
  */
  RefCountedPtr<BaseIF> m_plane;

  /*!
    @brief Noise function. 
  */
  RefCountedPtr<BaseIF> m_perlin;
};

#include <CD_NamespaceFooter.H>

#endif
