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

// Chombo includes
#include <REAL.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for encapsulating random number generation. 
  @details This class uses a singleton pattern to avoid code where many individual random number generators initialize RNGs using potentially the same seeds. This singleton
  is MPI-safe (ish). The user can specify a seed 's' where each MPI rank will initialize their RNG with seed 's + procID()'. Otherwise, a random seed is drawn by the master
  rank and broadcast to the other ranks. The rank seeds are like the above: 's+procID()'. 
*/
class Random {
public:

  /*!
    @brief Set the RNG seed.
    @details If using MPI, the seed is set to a_seed + procID().
    @param[in] a_seed Seed to use for RNG.
  */
  inline
  static void setSeed(const int a_seed);

  /*!
    @brief Set a random RNG seed
  */
  inline
  static void init();

  /*!
    @brief Get a uniform real number on the interval [0,1]
  */
  inline
  static Real getUniformReal01();

  /*!
    @brief Get a uniform real number on the interval [-1,1]
  */
  inline
  static Real getUniformReal11();

  /*!
    @brief Get a number from a normal distribution centered on zero and variance 1
  */
  inline
  static Real getNormal01();  

  /*!
    @brief Get a random direction in space.
    @details Uses Marsaglia algorithm. 
  */
  inline
  static RealVect getDirection();

  /*!
    @brief For getting a random number from a user-supplied distribution. T must be a distribution for which we can call T(s_rng)
    @param[in] a_distribution Distribution. Must have object of type
  */
  template<typename T>
  inline  
  static Real get( T& a_distribution);

protected:

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
