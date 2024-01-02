/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Random.H
  @brief  File containing some useful static methods related to random number generation
  @author Robert Marskar
*/

#ifndef CD_Random_H
#define CD_Random_H

// Std includes
#include <random>
#include <memory>
#include <mutex>
#include <type_traits>

// Chombo includes
#include <REAL.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for encapsulating random number generation. 
  @details The user can specify a seed 's' where each MPI rank will initialize their RNG with seed 's + procID()'. Note that unless the user specifies 
  Random.seed = <number> in the input script, the RNG will use a seed of '0 + procID()' for the various ranks. If the user specifies a <number> less than 0, 
  a random seed will be used.
*/
class Random
{
public:
  /*!
    @brief Disallowed constructor. 
  */
  Random() = delete;

  /*!
    @brief Deleted copy constructor. Must be impossible because this is a singleton. 
    @param[in] a_other Other RNG engine. 
  */
  Random(const Random& a_other) = delete;

  /*!
    @brief Deleted copy constructor. Must be impossible because this is a singleton. 
    @param[in] a_other Other RNG engine. 
  */
  Random(const Random&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment. 
    @param[in] a_other Other RNG engine. 
  */
  Random&
  operator=(const Random& a_other) = delete;

  /*!
    @brief Disallowed move assignment. 
    @param[in] a_other Other RNG engine. 
  */
  Random&
  operator=(const Random&& a_other) = delete;

  /*!
    @brief Get Poisson distributed number
    @param[in] a_mean Poisson mean value (inverse rate)
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  inline static T
  getPoisson(const Real a_mean);

  /*!
    @brief Get Poisson distributed number. 
    @param[in] a_N Number of trials
    @param[in] a_p The usual success probability in binomial distributions
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  inline static T
  getBinomial(const T a_N, const Real a_p) noexcept;

  /*!
    @brief Get a uniform real number on the interval [0,1]
  */
  inline static Real
  getUniformReal01();

  /*!
    @brief Get a uniform real number on the interval [-1,1]
  */
  inline static Real
  getUniformReal11();

  /*!
    @brief Get a number from a normal distribution centered on zero and variance 1
  */
  inline static Real
  getNormal01();

  /*!
    @brief Get a random direction in space.
    @details Uses Marsaglia algorithm. 
  */
  inline static RealVect
  getDirection();

  /*!
    @brief For getting a random number from a user-supplied distribution. T must be a distribution for which we can call T(s_rng)
    @param[in] a_distribution Distribution. Must have object of type
  */
  template <typename T>
  inline static Real
  get(T& a_distribution);

  /*!
    @brief For getting a random number from a user-supplied distribution. T must be a distribution for which we can call T(s_rng)
    @param[in] a_distribution Distribution. Must have object of type
  */
  template <typename T>
  inline static size_t
  getDiscrete(T& a_distribution);

  /*!
    @brief Seed the RNG
  */
  inline static void
  seed();

  /*!
    @brief Set the RNG seed.
    @details If using MPI, the seed is set to a_seed + procID().
    @param[in] a_seed Seed to use for RNG.
  */
  inline static void
  setSeed(const int a_seed);

  /*!
    @brief Set a random RNG seed
  */
  inline static void
  setRandomSeed();

  /*!
    @brief Return a random position in the cube (a_lo, a_hi);
    @param[in] a_lo Lower-left corner 
    @param[in] a_hi Upper-right corner 
  */
  inline static RealVect
  randomPosition(const RealVect a_lo, const RealVect a_hi) noexcept;

  /*!
    @brief Draw a random position physical somewhere in a grid cell.
    @details This is the version that respects cut-cells -- it will draw a position (through rejection sampling) such that the particle ends up 
    on the correct side of the EB. 
    @param[in] a_cellPos       Cell-center position.
    @param[in] a_lo            Lower-left corner of volume that encloses the cut-cell. This is relative to the unit cell. 
    @param[in] a_hi            Upper-right corner of volume that encloses the cut-cell. This is relative to the unit cell. 
    @param[in] a_bndryCentroid EB centroid position
    @param[in] a_normal        EB normal (points into the fluid region)
    @param[in] a_dx            Grid resolution
    @param[in] a_kappa         Cell volume
  */
  inline static RealVect
  randomPosition(const RealVect a_cellPos,
                 const RealVect a_lo,
                 const RealVect a_hi,
                 const RealVect a_bndryCentroid,
                 const RealVect a_normal,
                 const Real     a_dx,
                 const Real     a_kappa) noexcept;

  /*!
    @brief Draw a random position somewhere in a cut-cell
    @param[in] a_lo            Lower-left corner of volume that encloses the cut-cell. This is relative to the unit cell. 
    @param[in] a_hi            Upper-right corner of volume that encloses the cut-cell. This is relative to the unit cell. 
    @param[in] a_bndryCentroid EB centroid position
    @param[in] a_normal        EB normal (points into the fluid region)
    @return Returns a position inside the cut-cell and on the correct side of the EB. 
  */
  inline static RealVect
  randomPosition(const RealVect a_lo,
                 const RealVect a_hi,
                 const RealVect a_bndryCentroid,
                 const RealVect a_normal) noexcept;

private:
  /*!
    @brief For checking if RNG has been seeded or not.
  */
  static bool s_seeded;

  /*!
    @brief Random number generator
  */
  static std::mt19937_64 s_rng;

  /*!
    @brief For drawing random number on the interval [0,1]
  */
  static std::uniform_real_distribution<Real> s_uniform01;

  /*!
    @brief For drawing random number on the interval [-1,1]
  */
  static std::uniform_real_distribution<Real> s_uniform11;

  /*!
    @brief Normal distribution centered at zero with standard deviation of one
  */
  static std::normal_distribution<Real> s_normal01;
};

#include <CD_NamespaceFooter.H>

#include <CD_RandomImplem.H>

#endif
