/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSdf.H
  @brief  Declaration of signed distance Perlin noise
  @author Robert Marskar
*/

#ifndef CD_PerlinSdf_H
#define CD_PerlinSdf_H

// Std includes
#include <algorithm>

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that implements the improved Perlin noise function
  @details Typically, you will use this function to displace a level-set function by using SumIF. 
*/
class PerlinSdf : public BaseIF {
public:
  /*!
    @brief Full constructor
    @param[in] a_noiseAmp Noise amplitude
    @param[in] a_noiseFreq Base frequency
    @param[in] a_octaves Number of octaves to use for Perlin noise
    @param[in] a_persistence Reduction factor for each octave
    @param[in] a_reseed Reseed permutation vector (true) or use Ken Perlin's original one (false)
  */
  PerlinSdf(const Real     a_noiseAmp,
	    const RealVect a_noiseFreq,
	    const Real     a_persistence,
	    const int      a_octaves,
	    const bool     a_reseed = false);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input PerlinSdf
  */
  PerlinSdf(const PerlinSdf& a_inputIF);

  /*!
    @brief Random number generator
  */
  static int random(const int i);

  /*!
    @brief Destructor
  */
  virtual ~PerlinSdf();

  /*!
    @brief Level-set function
    @param[in] a_pos Position
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF* newImplicitFunction() const;
  
protected:

  /*!
    @brief Ken Perlin's original permutation arary
  */
  static int s_perlin[];

  /*!
    @brief Perlin noise frequency
  */
  RealVect m_noiseFreq;
  
  /*!
    @brief Perlin noise amplitude
  */
  Real m_noiseAmp;

  /*!
    @brief Perlin octave fade factor
  */
  Real m_persistence;

  /*!
    @brief Number of octaves to use for noise
  */
  int m_octaves;

  /*!
    @brief Permutation table. 
  */
  double p[512];

  /*!
    @brief Reseed function
  */
  virtual void reseed();
  
  /*!
    @brief Ken Perlin's 3D noise function. 
    @details For 2D, we use Ken Perlin's function with z = 0
    @param[in] a_x Position
    @param[in] a_y Position
    @param[in] a_z Position
  */
  double noise(const double a_x, const double a_y, const double a_z) const;

  /*!
    @brief Chombo interface to Perlin noise function
  */
  Real noise(const RealVect& a_pos) const;
  
  /*!
    @brief Octave noise function
  */
  Real octaveNoise(const RealVect& a_pos) const;
  
  /*!
    @brief Interpolation function
  */
  Real lerp(const Real t, const Real a, const Real b) const;

  /*!
    @brief Fade function
  */
  Real fade(const Real t) const;

  /*!
    @brief Gradient function using hash table lookup
  */
  Real grad(const int hash, const double x, const double y, const double z) const;
};

#include <CD_NamespaceFooter.H>


#endif
