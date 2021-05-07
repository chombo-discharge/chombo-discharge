/*!
  @file   MultiIndex.H
  @brief  Declaration of a SpaceDim multi-index type for use with higher order Taylor series
  @author Robert Marskar
  @date   Aug. 2017
*/

#ifndef _MULTIINDEX_
#define _MULTIINDEX_

#include <map>

#include "IntVect.H"
#include "RealVect.H"

/*!
  @brief SpaceDim multi-index type for use with higher order Taylor series
  @note The next() and hasNext() functions are just temporary functions; I'm just too lazy to write an iterator class (yet)
*/
class MultiIndex {
public:
  /*!
    @brief Weak construction is bad. 
  */
  MultiIndex() = delete;

  /*!
    @brief Copy constructor
  */
  MultiIndex(const MultiIndex& a_index) = default;

  /*!
    @brief Full constructor
    @param[i] a_order
  */
  MultiIndex(const int a_order);

  /*!
    @brief Destructor
  */
  ~MultiIndex() = default;

  /*!
    @brief Define function
  */
  void define(const int a_order);

  /*!
    @brief Factorial function
  */
  int factorial() const;

  /*!
    @brief Get index
  */
  IntVect getCurrentIndex() const;

  /*!
    @brief Get highest order for this index
  */
  int getOrder() const;

  /*!
    @brief Get the linear index (i.e. column number) for a given multi-index. 
  */
  int getLinearIndex(const IntVect a_MultiIndex) const;

  /*!
    @brief Get a multi-index corresponding to the linear index (i.e. column number)
  */
  IntVect getMultiIndex(const int a_linearIndex) const;

  /*!
    @brief Assignment operator
  */
  MultiIndex& operator=(const MultiIndex& a_other) = default;

  /*!
    @brief Get one index
  */
  int operator[](const int a_dir) const;

  /*!
    @brief Increment multi-index lexigraphically
  */
  void operator++() ;

  /*!
    @brief Reset multi-index to point to the first index
  */
  void reset();

  /*!
    @brief Check that multi-index is ok
  */
  bool ok() const;

  /*!
    @brief Get norm of multi-index.
    @details Returns |alpha| + |alpha2| + ....
  */
  int norm() const;

  /*!
    @brief Return number of indices in MultiIndex range
  */
  int getNumIndices() const;

  /*!
    @brief Power function for multi-indices
  */
  Real pow(const RealVect& a_vec);
  
protected:

  /*!
    @brief Maximum order representable by this multi-index
  */
  int m_order;

  /*!
    @brief Current idnex
  */
  IntVect m_curIndex;

  /*!
    @brief Pointer to index
  */
  std::vector<IntVect>::iterator m_iter;

  /*!
    @brief All indices
  */
  std::vector<IntVect> m_indices;

  /*!
    @brief Map for going from lexigraphically ordered multi-index to linear index
  */
  std::map<IntVect, int> m_mapToLinearIndex;

  /*!
    @brief Map going from linear index to multi-index
  */
  std::map<int, IntVect> m_mapToMultiIndex;

  /*!
    @brief Base factorial function
  */
  int factorial(const int a_n) const;

  /*!
    @brief Fill m_indices
  */
  void makeIndices();

  /*!
    @brief Fill maps
  */
  void makeMaps();

private:
  
  /*!
    @brief Get norm
  */
  int norm(const IntVect a_index) const;
};

#endif
