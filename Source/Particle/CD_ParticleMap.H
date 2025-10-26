/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleOpsImplem.H
  @brief  Implementation of CD_ParticleOps.H
  @author Robert Marskar
*/

#ifndef CD_ParticleMap_H
#define CD_ParticleMap_H

// Std includes
#include <cstdint>
#include <utility>
#include <unordered_map>
#include <map>

// Our includes
#include <CD_NamespaceHeader.H>

// Compile-time sanity checks (optional but nice)
static_assert(sizeof(uint32_t) == 4, "uint32_t must be 32 bits");
static_assert(sizeof(uint64_t) == 8, "uint64_t must be 64 bits");
static_assert(sizeof(size_t) == 8, "size_t must be 64 bits");

/*!
  @brief Function for encoding a pair of integers onto a 64-bit integer.
  @param[in] a_level Level index
  @param[in] a_gridIndex Grid index
  @details Used for hash maps.
*/
static inline uint64_t
encodePair(uint32_t a_level, uint32_t a_gridIndex) noexcept
{
  return (std::uint64_t(a_level) << 32) | std::uint64_t(a_gridIndex);
}

/*!
  @brief Function for decoding a 64-bit integer into two integers.
  @details Used for reversing hash maps.
  @param[in] key Hash key as 64-bit integer. 
*/
static inline std::pair<uint32_t, uint32_t>
decodePair(uint64_t key) noexcept
{
  return {uint32_t(key >> 32), uint32_t(key)};
}

/*!
  @brief Simple pair-hashing struct for interfacing into unordered_map
*/
struct PairHash
{
  std::size_t
  operator()(const std::pair<uint32_t, uint32_t>& p) const noexcept
  {
    return encodePair(p.first, p.second);
  }
};

/*!
  @brief Underlying particle map when gathering/scattering particles
*/
template <typename T>
#if 0
using ParticleMap = std::unordered_map<std::pair<uint32_t, uint32_t>, T, PairHash>;
#else
using ParticleMap = std::map<std::pair<unsigned int, unsigned int>, T>;
#endif

#include <CD_NamespaceFooter.H>

#endif
