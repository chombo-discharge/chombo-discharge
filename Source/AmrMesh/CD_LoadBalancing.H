/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_LoadBalancing.H
  @brief  Declaration of a static class for various load balancing operations.
  @author Robert Marskar
*/

#ifndef CD_LOADBALANCING_H
#define CD_LOADBALANCING_H

// Std includes
#include <cstdint>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_BoxSorting.H>
#include <CD_Loads.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Static class which contains various load balancing routines. 
*/
class LoadBalancing
{
public:
  /**
    @brief Load balancing, assigning ranks to boxes. 
    @details This version does not take into loads already assigned to ranks and will thus always start assigning from rank 0. 
    @param[out] a_ranks Vector containing processor IDs corresponding to boxes (and loads)
    @param[in]  a_loads Computational loads
    @param[in]  a_boxes Grid boxes
  */
  template <class T>
  static void
  makeBalance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes);

  /**
    @brief Load balancing, assigning ranks to boxes. 
    @details This version takes into account loads already assigned to ranks, and assigns grid subsets to the least
    loaded ranks. 
    @param[out]   a_ranks Vector containing processor IDs corresponding to boxes (and loads)
    @param[in,out] a_rankLoads MPI rank loads so far
    @param[in]    a_boxLoads Computational loads for each box
    @param[in]    a_boxes Grid boxes
  */
  template <class T>
  static void
  makeBalance(Vector<int>& a_ranks, Loads& a_rankLoads, const Vector<T>& a_boxLoads, const Vector<Box>& a_boxes);

  /**
    @brief Sorts boxes and loads over a hierarchy according to some sorting criterion.
    @param[in,out] a_boxes Grid boxes
    @param[in,out] a_loads Computational loads
    @param[in]    a_whichSorting
  */
  template <typename T>
  static void
  sort(Vector<Vector<Box>>& a_boxes, Vector<Vector<T>>& a_loads, const BoxSorting a_whichSorting);

  /**
    @brief Sorting function for sorting boxes and loads
    @param[in,out] a_boxes Grid boxes
    @param[in,out] a_loads Computational loads
    @param[in] a_which Which
  */
  template <typename T>
  static void
  sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const BoxSorting a_which);

  /**
    @brief Sorts boxes (according to sorting argument)
    @param[in,out] a_boxes Grid boxes
    @param[in] a_which Description
  */
  static void
  sort(Vector<Box>& a_boxes, const BoxSorting a_which);

  /**
    @brief Gather boxes
    @param[in,out] a_boxes Local boxes
  */
  static void
  gatherBoxes(Vector<Box>& a_boxes);

  /**
    @brief Gather computational loads
    @param[in,out] a_loads Computational loads
  */
  static void
  gatherLoads(Vector<Real>& a_loads);

  /**
    @brief Gather compute loads
    @param[in,out] a_loads Computational loads
  */
  static void
  gatherLoads(Vector<int>& a_loads);

  /**
    @brief Gather compute loads
    @param[in,out] a_loads Computational loads
  */
  static void
  gatherLoads(Vector<long>& a_loads);

  /**
    @brief Gather boxes and loads
    @param[in,out] a_boxes Grid boxes
    @param[in,out] a_loads Computational loads
  */
  static void
  gatherBoxesAndLoads(Vector<Box>& a_boxes, Vector<int>& a_loads);

protected:
  /**
    @brief Utility function which packs boxes and loads into a vector of pairs
    @param[in] a_boxes Grid boxes
    @param[in] a_loads Computational loads
    @return Return value
  */
  template <class T>
  static std::vector<std::pair<Box, T>>
  packPairs(const Vector<Box>& a_boxes, const Vector<T>& a_loads);

  /**
    @brief Splits vector pair into separate boxes and loads. 
    @param[out] a_boxes Boxes.
    @param[out] a_loads Loads.
    @param[out] a_pairs Boxes and loads, stored as a pair. 
  */
  template <class T>
  static void
  unpackPairs(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T>>& a_pairs);

  /**
    @brief Standard box sorting, calls C++ std::sort
    @param[in,out] a_boxes Boxes to be sorted. 
    @param[in,out] a_loads Loads to be sorted. 
    @details This uses lexicographical indexing of boxes. On output, the boxes/loads are sorted with this ordering. 
  */
  template <class T>
  static void
  standardSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /**
    @brief Randomly shuffles boxes and loads.
    @param[in,out] a_boxes Grid boxes to be shuffled. 
    @param[in,out] a_loads Computational loads to be shuffled. 
    @details This shuffles the boxes and loads so that the box distribution is random. 
  */
  template <class T>
  static void
  shuffleSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /**
    @brief Generic SFC sorting function
    @param[in,out] a_boxes Grid boxes to be sorted
    @param[in,out] a_loads Computational loads to be sorted along with the boxes.
    @param[in]    a_sfcEncoder Function for turning a grid index into an SFC index.
  */
  template <class T, int DIM>
  static void
  sortSFC(Vector<Box>&                                                    a_boxes,
          Vector<T>&                                                      a_loads,
          const std::function<uint64_t(const std::array<uint32_t, DIM>)>& a_sfcEncoder) noexcept;

public:
  /**
    @brief Compute a Morton (Z-order) code from integer coordinates, using 21 bits per direction.
    @param[in] a_coords Integer coordinates (each < 2^21).
    @return The interleaved Morton index.
  */
  template <int DIM>
  static uint64_t
  mortonIndex(const std::array<uint32_t, DIM>& a_coords) noexcept;

  /**
    @brief Compute a Hilbert-curve index from integer coordinates, using 21 bits per direction.
    @param[in] a_coords Integer coordinates (each < 2^21).
    @return The Hilbert index.
  */
  template <int DIM>
  static uint64_t
  hilbertIndex(const std::array<uint32_t, DIM>& a_coords);
};

#include <CD_NamespaceFooter.H>

#include <CD_LoadBalancingImplem.H>

#endif
