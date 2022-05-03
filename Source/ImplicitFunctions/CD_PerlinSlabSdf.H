/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSlabSdf.H
  @brief  Declaration of a noisy slab, i.e. a box with surface roughness
  @author Robert Marskar
*/

#ifndef CD_PerlinSlabSdf_H
#define CD_PerlinSlabSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Implicit function class describing a slab (box) with surface noise on the top. 
*/
class PerlinSlabSdf : public BaseIF
{
public:
  /*!
    @brief No weak constructor allowed, use the full constructor
  */
  PerlinSlabSdf() = delete;

  /*!
    @brief Constructor
    @param[in] a_ccPoint     Center point on the slab side with surface noise
    @param[in] a_normal      Normal vector of slab side with surface noise
    @param[in] a_xyz         Width, depth, and height of slab
    @param[in] a_curv        Corner curvatures
    @param[in] a_noiseFreq   Noise frequency
    @param[in] a_octaves     Number of octaves
    @param[in] a_noiseAmp    Noise amplitude
    @param[in] a_persistence Factor for noise octaves. 
    @param[in] a_cornerCurv  Corner curvatures
    @param[in] a_reseed      If true, reseed the RNG. Otherwise, the the original Perlin hash table. 
    @note The persistence argument indicates how we add octaves. With 0.5 we use half the amplitude and twice the frequency for the next octave. 
    @param[in] a_fluidInside Fluid on inside or not
  */
  PerlinSlabSdf(const RealVect a_ccPoint,
                const RealVect a_normal,
                const RealVect a_xyz,
                const RealVect a_noiseFreq,
                const int      a_octaves,
                const Real     a_noiseAmp,
                const Real     a_persistence,
                const Real     a_cornerCurv,
                const bool     a_reseed,
                const bool     a_fluidInside);

  /*!
    @brief Copy constructor. Copies over implicit function.
  */
  PerlinSlabSdf(const PerlinSlabSdf& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~PerlinSlabSdf();

  /*!
    @brief Value function. Returns distance to object. 
    @param[in] a_pos Position. 
  */
  virtual Real
  value(const RealVect& a_pos) const override;

  /*!
    @brief Factory function
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief Implicit function. This is a union between cutting planes and a perlin plane
  */
  RefCountedPtr<BaseIF> m_baseif;

  /*!
    @brief Fluid on inside or not
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
