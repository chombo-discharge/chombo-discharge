/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LoadBalancing.H
  @brief  Declaration of a static class for various load balancing operations.
  @author Robert Marskar
*/

#ifndef CD_LoadBalancing_H
#define CD_LoadBalancing_H

// Std includes
#include <cstdint>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_BoxSorting.H>
#include <CD_Loads.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Static class which contains various load balancing routines. 
*/
class LoadBalancing
{
public:
  /*!
    @brief Load balancing, assigning ranks to boxes. 
    @details This version does not take into loads already assigned to ranks and will thus always start assigning from rank 0. 
    @param[out] a_ranks Vector containing processor IDs corresponding to boxes (and loads)
    @param[in]  a_loads Computational loads
    @param[in]  a_boxes Grid boxes
  */
  template <class T>
  static void
  makeBalance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes);

  /*!
    @brief Load balancing, assigning ranks to boxes. 
    @details This version takes into account loads already assigned to ranks, and assigns grid subsets to the least
    loaded ranks. 
    @param[out]   a_ranks Vector containing processor IDs corresponding to boxes (and loads)
    @param[inout] a_rankLoads MPI rank loads so far
    @param[in]    a_boxLoads Computational loads for each box
    @param[in]    a_boxes Grid boxes
  */
  template <class T>
  static void
  makeBalance(Vector<int>& a_ranks, Loads& a_rankLoads, const Vector<T>& a_boxLoads, const Vector<Box>& a_boxes);

  /*!
    @brief Sorts boxes and loads over a hierarchy according to some sorting criterion.
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @param[in]    a_whichSorting
    @return Returns a_boxes and a_loads, sorting according to the argument a_whichSorting.
  */
  template <typename T>
  static void
  sort(Vector<Vector<Box>>& a_boxes, Vector<Vector<T>>& a_loads, const BoxSorting a_whichSorting);

  /*!
    @brief Sorting function for sorting boxes and loads
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @param[in]    a_whichSorting
    @return Returns a_boxes and a_loads, sorting according to the argument a_whichSorting.
  */
  template <typename T>
  static void
  sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const BoxSorting a_which);

  /*!
    @brief Sorts boxes (according to sorting argument)
    @param[inout] a_boxes Grid boxes
  */
  static void
  sort(Vector<Box>& a_boxes, const BoxSorting a_which);

  /*!
    @brief Gather boxes
    @param[inout] a_boxes Local boxes
    @return On output, a_boxes contains the global boxes. 
  */
  static void
  gatherBoxes(Vector<Box>& a_boxes);

  /*!
    @brief Gather computational loads
    @param[inout] a_loads Computational loads
    @return On output, a_loads contains the loads from all ranks. 
  */
  static void
  gatherLoads(Vector<Real>& a_loads);

  /*!
    @brief Gather compute loads
    @param[inout] a_loads Computational loads
    @return On output, a_loads contains the loads from all ranks. 
  */
  static void
  gatherLoads(Vector<int>& a_loads);

  /*!
    @brief Gather compute loads
    @param[inout] a_loads Computational loads
    @return On output, a_loads contains the loads from all ranks. 
  */
  static void
  gatherLoads(Vector<long>& a_loads);

  /*!
    @brief Gather boxes and loads
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @return On output, a_boxes and a_loads have been globally gathered from all ranks. 
  */
  static void
  gatherBoxesAndLoads(Vector<Box>& a_boxes, Vector<int>& a_loads);

protected:
  /*!
    @brief Utility function which packs boxes and loads into a vector of pairs
    @param[in] a_boxes Grid boxes
    @param[in] a_loads Computational loads
  */
  template <class T>
  static std::vector<std::pair<Box, T>>
  packPairs(const Vector<Box>& a_boxes, const Vector<T>& a_loads);

  /*!
    @brief Splits vector pair into separate boxes and loads. 
    @param[out] a_boxes Boxes.
    @param[out] a_loads Loads.
    @param[out] a_pairs Boxes and loads, stored as a pair. 
  */
  template <class T>
  static void
  unpackPairs(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T>>& a_pairs);

  /*!
    @brief Standard box sorting, calls C++ std::sort
    @param[inout] a_boxes Boxes to be sorted. 
    @param[inout] a_loads Loads to be sorted. 
    @details This uses lexicographical indexing of boxes. On output, the boxes/loads are sorted with this ordering. 
  */
  template <class T>
  static void
  standardSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Randomly shuffles boxes and loads.
    @param[inout] a_boxes Grid boxes to be shuffled. 
    @param[inout] a_loads Computational loads to be shuffled. 
    @details This shuffles the boxes and loads so that the box distribution is random. 
  */
  template <class T>
  static void
  shuffleSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Sort boxes using Morton code.
    @param[inout] a_boxes Grid boxes to be sorted. 
    @param[inout] a_loads Computational loads to be sorted.
    @details On output, a_boxes and a_loads are Morton sorted. 
  */
  template <class T>
  static void
  mortonSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Compute Morton code using 21 bits per direction. 
    @param[in] a_p Integer index
  */
  static uint64_t
  mortonCode21(const IntVect a_iv) noexcept;

  /*!
    @brief Sort boxes using Hilbert code.
    @param[inout] a_boxes Grid boxes to be sorted. 
    @param[inout] a_loads Computational loads to be sorted.
    @details On output, a_boxes and a_loads are Hilbert sorted. 
  */
  template <class T>
  static void
  hilbertSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Helper function for computing a Hilbert index
  */
  template <int DIM>
  static uint64_t
  hilbertIndex(const std::array<uint32_t, DIM>& a_coords);
};

#include <CD_NamespaceFooter.H>

#include <CD_LoadBalancingImplem.H>

#endif
