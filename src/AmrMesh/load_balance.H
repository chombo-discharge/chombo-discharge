/*!
  @file load_balance.H
  @brief Declaration of various load-balancing routines
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _CD_LoadBalance_H_
#define _CD_LoadBalance_H_

#include "mfis.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Enum for sorting boxes
*/
enum class box_sorting {
  none,
  std,
  shuffle,
  morton
};

/*!
  @brief Static class that does various kinds of load-balancing
*/
class load_balance {
public:
  
  /*!
    @brief Load balancing
  */
  template <class T>
  static void make_balance(Vector<int>& a_ranks, const Vector<T>& a_loads, const Vector<Box>& a_boxes);

  /*!
    @brief Load balancing function which load balances level by level. 
  */
  template <class T>
  static void level_by_level(Vector<Vector<int> >&       a_procs,
			     const Vector<Vector<T> >&   a_loads,
			     const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Load balancing function which load balances all boxes at once, not caring about actual levels. 
  */
  template <class T>
  static void hierarchy(Vector<Vector<int> >&       a_procs,
			const Vector<Vector<T> >&   a_loads,
			const Vector<Vector<Box> >& a_boxes);

  /*!
    @brief Sort boxes and loads over a hierarchy
  */
  template <typename T>
  static void sort(Vector<Vector<Box> >& a_boxes, Vector<Vector<T> >& a_loads, const box_sorting a_which);
  
  /*!
    @brief Main sorting function for sorting boxes and loads
  */
  template <typename T>
  static void sort(Vector<Box>& a_boxes, Vector<T>& a_loads, const box_sorting a_which);

  /*!
    @brief Load balancing. Loads are box volumes. 
  */
  static void make_balance(Vector<int>& a_levelRanks, const Vector<Box>& a_levelBoxes);

  /*!
    @brief Round robin distribution of boxes
  */
  static void round_robin(Vector<int>& a_ranks, const Vector<Box>& a_boxes);

  /*!
    @brief Main sorting function for sorting boxes and loads
    @details This just calls the other version. 
  */
  static void sort(Vector<Box>& a_boxes, const box_sorting a_which);

  /*!
    @brief Gather boxes
  */
  static void gather_boxes(Vector<Box>& a_boxes);

  /*!
    @brief Gather compute loads
  */
  static void gather_loads(Vector<Real>& a_loads);

  /*!
    @brief Gather compute loads
  */
  static void gather_loads(Vector<int>& a_loads);

  /*!
    @brief Gather boxes and loads
  */
  static void gather_boxes_and_loads(Vector<Box>& a_boxes, Vector<int>& a_loads);

protected:

  /*!
    @brief Make a vector pair
  */
  template <class T>
  static std::vector<std::pair<Box, T> > packPairs(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Splits vector pair
  */
  template <class T>
  static void unpackPairs(Vector<Box>& a_boxes, Vector<T>& a_loads, const std::vector<std::pair<Box, T> >& a_pairs);

  /*!
    @brief Standard box sorting, calls C++ std::sort
  */
  template <class T>
  static void std_sort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Random shuffle boxes
  */
  template <class T>
  static void shuffle_sort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Morton sort boxes
  */
  template <class T>
  static void morton_sort(Vector<Box>& a_boxes, Vector<T>& a_loads);

  /*!
    @brief Morton comparator
  */
  template <class T>
  static bool mortonComparator(const int a_maxBits, const std::pair<Box,T>& a_lhs, const std::pair<Box, T>& a_rhs);

  /*!
    @brief Get maximum bits
  */
  static int maxBits(std::vector<Box>::iterator a_first, std::vector<Box>::iterator a_last);

};
#include "CD_NamespaceFooter.H"

#include "load_balanceI.H"

#endif
