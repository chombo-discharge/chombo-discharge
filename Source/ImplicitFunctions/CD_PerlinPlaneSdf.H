/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinPlaneSdf.H
  @brief  Declaration of a signed distance function for a noisy plane
  @author Robert Marskar
*/

#ifndef CD_PerlinPlaneSdf_H
#define CD_PerlinPlaneSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Infinite plane with surface noise (Perlin noise). 
*/
class PerlinPlaneSdf : public BaseIF
{
public:
  /*!
    @brief Full constructor
    @param[in] a_point       Point on the plane
    @param[in] a_normal      Plane normal vector
    @param[in] a_inside      Inside or outside
    @param[in] a_noiseAmp    Noise amplitude
    @param[in] a_noiseFreq   Noise frequency
    @param[in] a_persistence Factor for noise octaves. 
    @param[in] a_octaves     Number of octaves
    @param[in] a_reseed      If true, reseed the RNG. Otherwise, the the original Perlin hash table. 
    @note The persistence argument indicates how we add octaves. With 0.5 we use half the amplitude and twice the frequency for the next octave. 
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
    @param[in] a_inputIF Other plane
  */
  PerlinPlaneSdf(const PerlinPlaneSdf& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~PerlinPlaneSdf();

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
    @brief Plane normal
  */
  RealVect m_normal;

  /*!
    @brief Plane point
  */
  RealVect m_point;

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
