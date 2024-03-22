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
  @details Typically, you will use this function to displace a level-set function by using SumIF. However, Perlin
  noise is also a signed distance function, and so it can be used as an implicit function as well. 
  @note See the original paper by Ken Perlin for understanding the algorithm: "Improving Noise. Ken Perlin (2002)"
*/
class PerlinSdf : public BaseIF
{
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
    @brief Simple random number generator
    @param[in] seed
  */
  static int
  random(const int i);

  /*!
    @brief Destructor
  */
  virtual ~PerlinSdf();

  /*!
    @brief Level-set function
    @param[in] a_pos Position
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /*!
    @brief Factory method
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /*!
    @brief Ken Perlin's original permutation array
  */
  constexpr static int m_permutationTable[256] = {
    151, 160, 137, 91,  90,  15,  131, 13,  201, 95,  96,  53,  194, 233, 7,   225, 140, 36,  103, 30,  69,  142,
    8,   99,  37,  240, 21,  10,  23,  190, 6,   148, 247, 120, 234, 75,  0,   26,  197, 62,  94,  252, 219, 203,
    117, 35,  11,  32,  57,  177, 33,  88,  237, 149, 56,  87,  174, 20,  125, 136, 171, 168, 68,  175, 74,  165,
    71,  134, 139, 48,  27,  166, 77,  146, 158, 231, 83,  111, 229, 122, 60,  211, 133, 230, 220, 105, 92,  41,
    55,  46,  245, 40,  244, 102, 143, 54,  65,  25,  63,  161, 1,   216, 80,  73,  209, 76,  132, 187, 208, 89,
    18,  169, 200, 196, 135, 130, 116, 188, 159, 86,  164, 100, 109, 198, 173, 186, 3,   64,  52,  217, 226, 250,
    124, 123, 5,   202, 38,  147, 118, 126, 255, 82,  85,  212, 207, 206, 59,  227, 47,  16,  58,  17,  182, 189,
    28,  42,  223, 183, 170, 213, 119, 248, 152, 2,   44,  154, 163, 70,  221, 153, 101, 155, 167, 43,  172, 9,
    129, 22,  39,  253, 19,  98,  108, 110, 79,  113, 224, 232, 178, 185, 112, 104, 218, 246, 97,  228, 251, 34,
    242, 193, 238, 210, 144, 12,  191, 179, 162, 241, 81,  51,  145, 235, 249, 14,  239, 107, 49,  192, 214, 31,
    181, 199, 106, 157, 184, 84,  204, 176, 115, 121, 50,  45,  127, 4,   150, 254, 138, 236, 205, 93,  222, 114,
    67,  29,  24,  72,  243, 141, 128, 195, 78,  66,  215, 61,  156, 180};

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
  virtual void
  reseed();

  /*!
    @brief Ken Perlin's 3D noise function. 
    @details For 2D, we use Ken Perlin's function with z = 0
    @param[in] a_x Position
    @param[in] a_y Position
    @param[in] a_z Position
  */
  double
  noise(const double a_x, const double a_y, const double a_z) const;

  /*!
    @brief Chombo interface to Perlin noise function
  */
  Real
  noise(const RealVect& a_pos) const;

  /*!
    @brief Octave noise function. Returns noise over over several octaves. 
  */
  Real
  octaveNoise(const RealVect& a_pos) const;

  /*!
    @brief Interpolation function
  */
  Real
  lerp(const Real t, const Real a, const Real b) const;

  /*!
    @brief Fade function
  */
  Real
  fade(const Real t) const;

  /*!
    @brief Gradient function using hash table lookup
  */
  Real
  grad(const int hash, const double x, const double y, const double z) const;
};

#include <CD_NamespaceFooter.H>

#endif
