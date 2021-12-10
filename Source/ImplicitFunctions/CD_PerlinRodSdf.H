/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  @brief A rod geometry, but with noise on one of the end caps.
*/
class PerlinRodSdf : public PerlinSphereSdf {
public:

  /*!
    @brief No weak construction
  */
  PerlinRodSdf() = delete;
  
  /*!
    @brief Constructor
    @param[in] a_rad         Rod radius
    @param[in] a_center1     One end of the rod
    @param[in] a_center2     Other end of the rod. 
    @param[in] a_inside      Inside or outside
    @param[in] a_noiseAmp    Noise amplitude
    @param[in] a_noiseFreq   Noise frequency
    @param[in] a_persistence Factor for noise octaves. 
    @param[in] a_octaves     Number of octaves
    @param[in] a_reseed      If true, reseed the RNG. Otherwise, the the original Perlin hash table. 
    @note The persistence argument indicates how we add octaves. With 0.5 we use half the amplitude and twice the frequency for the next octave. 
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
    @param[in] a_inputIF Other object
  */
  PerlinRodSdf(const PerlinRodSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~PerlinRodSdf();

  /*!
    @brief Value function
    @param[in] a_pos Position
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
