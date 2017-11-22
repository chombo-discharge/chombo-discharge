/*!
  @file perlin_if.H
  @brief Declaration of a perlin_if function
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PERLIN_IF_
#define _PERLIN_IF_

#include <BaseIF.H>
#include <algorithm>

/*!
  @class perlin_if
  @brief Class that implements an improved Perlin noise function
  @details Typically, you will use this function to displace a level-set function by using SumIF. 
*/
class perlin_if : public BaseIF {
public:
  /*!
    @brief Full constructor
    @param[in] a_noiseAmp Noise amplitude
    @param[in] a_noiseFreq Base frequency
    @param[in] a_octaves Number of octaves to use for Perlin noise
    @param[in] a_persistence Reduction factor for each octave
    @param[in] a_reseed Reseed permutation vector (true) or use Ken Perlin's original one (false)
  */
  perlin_if(const Real&     a_noiseAmp,
	    const RealVect& a_noiseFreq,
	    const Real&     a_persistence,
	    const int&      a_octaves,
	    const bool&     a_reseed = false);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input perlin_if
  */
  perlin_if(const perlin_if& a_inputIF);

  /*!
    @brief Random number generator
  */
  static int random(const int i);

  /*!
    @brief Destructor
  */
  virtual ~perlin_if();

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
  const double noise(const double a_x, const double a_y, const double a_z) const;

  /*!
    @brief Chombo interface to Perlin noise function
  */
  const Real noise(const RealVect& a_pos) const;
  
  /*!
    @brief Octave noise function
  */
  const Real octaveNoise(const RealVect& a_pos) const;
  
  /*!
    @brief Interpolation function
  */
  const Real lerp(const Real t, const Real a, const Real b) const;

  /*!
    @brief Fade function
  */
  const Real fade(const Real t) const;

  /*!
    @brief Gradient function using hash table lookup
  */
  const Real grad(const int hash, const double x, const double y, const double z) const;
};


#endif
