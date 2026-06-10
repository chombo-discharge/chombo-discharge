/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MultiIndex.H
  @brief  Declaration of a SpaceDim multi-index type for use with higher order Taylor series
  @author Robert Marskar
*/

#ifndef CD_MULTIINDEX_H
#define CD_MULTIINDEX_H

// Std includes
#include <map>

// Chombo includes
#include <IntVect.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief SpaceDim multi-index type for use with higher order Taylor series
  @note The next() and hasNext() functions are just temporary functions; I'm just too lazy to write an iterator class (yet)
*/
class MultiIndex
{
public:
  /**
    @brief Weak construction is bad. 
  */
  MultiIndex() = delete;

  /**
    @brief Copy constructor
    @param[in,out] a_index Description
  */
  MultiIndex(const MultiIndex& a_index) = default;

  /**
    @brief Full constructor
    @param[in] a_order Multi-index order
  */
  MultiIndex(const int a_order);

  /**
    @brief Destructor
  */
  ~MultiIndex() = default;

  /**
    @brief Define function
    @param[in] a_order Description
  */
  void
  define(const int a_order);

  /**
    @brief Factorial function
    @return Return value
  */
  int
  factorial() const;

  /**
    @brief Get index
    @return Current index
  */
  IntVect
  getCurrentIndex() const;

  /**
    @brief Get highest order for this index
    @return Order
  */
  int
  getOrder() const;

  /**
    @brief Get the linear index (i.e. column number) for a given multi-index. 
    @return Linear index
    @param[in] a_MultiIndex Description
  */
  int
  getLinearIndex(const IntVect a_MultiIndex) const;

  /**
    @brief Get a multi-index corresponding to the linear index (i.e. column number)
    @return Multi index
    @param[in] a_linearIndex Description
  */
  IntVect
  getMultiIndex(const int a_linearIndex) const;

  /**
    @brief Assignment operator
    @return Return value
    @param[in,out] a_other Description
  */
  MultiIndex&
  operator=(const MultiIndex& a_other) = default;

  /**
    @brief Get one index
    @return Return value
    @param[in] a_dir Description
  */
  int
  operator[](const int a_dir) const;

  /**
    @brief Increment multi-index lexigraphically
  */
  void
  operator++();

  /**
    @brief Reset multi-index to point to the first index
  */
  void
  reset();

  /**
    @brief Check that multi-index is ok
    @return Return value
  */
  bool
  ok() const;

  /**
    @brief Get norm of multi-index.
    @details Returns |alpha| + |alpha2| + ....
    @return Return value
  */
  int
  norm() const;

  /**
    @brief Return number of indices in MultiIndex range
    @return Num indices
  */
  int
  getNumIndices() const;

  /**
    @brief Power function for multi-indices
    @return Return value
    @param[in,out] a_vec Description
  */
  Real
  pow(const RealVect& a_vec);

protected:
  /**
    @brief Maximum order representable by this multi-index
  */
  int m_order;

  /**
    @brief Current idnex
  */
  IntVect m_curIndex;

  /**
    @brief Pointer to index
  */
  std::vector<IntVect>::iterator m_iter;

  /**
    @brief All indices
  */
  std::vector<IntVect> m_indices;

  /**
    @brief Map for going from lexigraphically ordered multi-index to linear index
  */
  std::map<IntVect, int> m_mapToLinearIndex;

  /**
    @brief Map going from linear index to multi-index
  */
  std::map<int, IntVect> m_mapToMultiIndex;

  /**
    @brief Base factorial function
    @return Return value
    @param[in] a_n Description
  */
  int
  factorial(const int a_n) const;

  /**
    @brief Fill m_indices
  */
  void
  makeIndices();

  /**
    @brief Fill maps
  */
  void
  makeMaps();

private:
  /**
    @brief Get norm
  */
  static int
  norm(const IntVect a_iv);
};

#include <CD_NamespaceFooter.H>

#endif
