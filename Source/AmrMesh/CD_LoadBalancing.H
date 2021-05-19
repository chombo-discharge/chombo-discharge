/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   LoadBalancing.H
  @brief  Declaration of a static class for various load balancing operations.
  @author Robert Marskar
*/

#ifndef CD_LoadBalancing_H
#define CD_LoadBalancing_H

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_BoxSorting.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace containing various load balancing routines.
*/
class LoadBalancing {
public:
  
  /*!
    @brief Load balancing assigning ranks to boxes. 
    @param[out] a_ranks Vector containing processor IDs corresponding to boxes (and loads)
    @param[in]  a_loads Computational loads
    @param[in]  a_boxes Grid boxes
  */
  template <class T>
  static void makeBalance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes);

  /*!
    @brief Load balancing function which load balances level by level. 
    @param[out] a_procs Vector containing processor IDs corresponding to boxes (and loads)
    @param[in]  a_loads Computational loads
    @param[in]  a_boxes Grid boxes
  */
  template <class T>
  static void balanceLevelByLevel(Vector<Vector<int> >&       a_procs,
				  const Vector<Vector<T> >&   a_loads,
				  const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Load balancing function which load balances all boxes at once, not caring about actual levels. 
    @param[out] a_procs Vector containing processor IDs corresponding to boxes (and loads)
    @param[in]  a_loads Computational loads
    @param[in]  a_boxes Grid boxes
  */
  template <class T>
  static void hierarchy(Vector<Vector<int> >&       a_procs,
			const Vector<Vector<T> >&   a_loads,
			const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Sorts boxes and loads over a hierarchy according to some sorting criterion.
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @param[in]    a_whichSorting
    @return Returns a_boxes and a_loads, sorting according to the argument a_whichSorting.
  */
  template <typename T>
  static void sort(Vector<Vector<Box> >& a_boxes, Vector<Vector<T> >& a_loads, const BoxSorting a_whichSorting);
  
  /*!
    @brief Sorting function for sorting boxes and loads
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @param[in]    a_whichSorting
    @return Returns a_boxes and a_loads, sorting according to the argument a_whichSorting.
  */
  template <typename T>
  static void sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const BoxSorting a_which);

  /*!
    @brief Load balancing of grid boxes using the box volume as computational load
    @param[out] a_ranks Processor IDs
    @param[in]  a_boxes Grid boxes
  */
  static void makeBalance(Vector<int>& a_ranks, const Vector<Box>& a_boxes);

  /*!
    @brief Round robin distribution of boxes
    @param[out] a_ranks Processor IDs
    @param[in]  a_boxes Grid boxes
  */
  static void roundRobin(Vector<int>& a_ranks, const Vector<Box>& a_boxes);

  /*!
    @brief Sorts boxes (according to sorting argument)
    @param[inout] a_boxes Grid boxes
  */
  static void sort(Vector<Box>& a_boxes, const BoxSorting a_which);

  /*!
    @brief Gather boxes
    @param[inout] a_boxes Local boxes
    @return On output, a_boxes contains the global boxes. 
  */
  static void gatherBoxes(Vector<Box>& a_boxes);

  /*!
    @brief Gather computational loads
    @param[inout] a_loads Computational loads
    @return On output, a_loads contains the loads from all ranks. 
  */
  static void gatherLoads(Vector<Real>& a_loads);

  /*!
    @brief Gather compute loads
    @param[inout] a_loads Computational loads
    @return On output, a_loads contains the loads from all ranks. 
  */
  static void gatherLoads(Vector<int>& a_loads);

  /*!
    @brief Gather boxes and loads
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
    @return On output, a_boxes and a_loads have been globally gathered from all ranks. 
  */
  static void gatherBoxes_and_loads(Vector<Box>& a_boxes, Vector<int>& a_loads);

protected:

  /*!
    @brief Utility function which packs boxes and loads into a vector of pairs
    @param[in] a_boxes Grid boxes
    @param[in] a_loads Computational loads
  */
  template <class T>
  static std::vector<std::pair<Box, T> > packPairs(const Vector<Box>& a_boxes, const Vector<T>& a_loads);

  /*!
    @brief Splits vector pair
  */
  template <class T>
  static void unpackPairs(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T> >& a_pairs);

  /*!
    @brief Standard box sorting, calls C++ std::sort
  */
  template <class T>
  static void standardSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Randomly shuffles boxes and loads.
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
  */
  template <class T>
  static void shuffleSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Sort boxes using Morton code.
    @param[inout] a_boxes Grid boxes
    @param[inout] a_loads Computational loads
  */
  template <class T>
  static void mortonSort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Morton comparator
    @param[in] a_maxBits Maximum bits
    @param[in] a_lhs     Left hand side
    @param[in] a_rhs     Right hand side
    @return Returns true if a_lhs has a lower Morton index than a_rhs
  */
  template <class T>
  static bool mortonComparator(const int a_maxBits, const std::pair<Box,T>& a_lhs, const std::pair<Box, T>& a_rhs);

  /*!
    @brief Utility function for Morton sorting which needs the maximum bits
    @param[in] a_first Iterator 
    @param[in] a_second Iterator 
  */
  static int maxBits(std::vector<Box>::iterator a_first, std::vector<Box>::iterator a_last);
};

#include <CD_NamespaceFooter.H>

#include <CD_LoadBalancingImplem.H>

#endif
