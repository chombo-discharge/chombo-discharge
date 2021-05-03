/*!
  @file   MultiIndex.H
  @brief  Declaration of a SpaceDim multi-index type for use with higher order Taylor series
  @author Robert Marskar
  @date   Aug. 2017
*/

#ifndef _MULTIINDEX_
#define _MULTIINDEX_

#include "IntVect.H"
#include "RealVect.H"

/*!
  @brief SpaceDim multi-index type for use with higher order Taylor series
  @note The next() and hasNext() functions are just temporary functions; I'm just too lazy to write an iterator class (yet)
*/
class MultiIndex {
public:
  /*!
    @brief Weak construction
  */
  MultiIndex();

  /*!
    @brief Full define function
  */
  MultiIndex(const IntVect& a_index);

  /*!
    @brief Copy constructor
  */
  MultiIndex(const MultiIndex& a_index);

  /*!
    @brief Destructor
  */
  ~MultiIndex();

  /*!
    @brief Define function
    @param[in] a_index Index
  */
  void define(const IntVect& a_index);

  /*!
    @brief Factorial function
  */
  int factorial() const;

  /*!
    @brief Assigment operator
  */
  MultiIndex& operator=(const MultiIndex& a_index);

  /*!
    @brief Get index
  */
  IntVect getIntVect() const;

  /*!
    @brief Get one index
  */
  int operator[](const int a_dir) const;

  /*!
    @brief Equal operator 
  */
  bool operator==(const MultiIndex& a_index) const;

  /*!
    @brief Equal operator 
  */
  bool operator==(const int a_Q) const;

  /*!
    @brief Comparison operator. This uses lexicographical ordering
  */
  bool operator>(const MultiIndex& a_index) const;
    /*!
    @brief Comparison operator. Returns true if abs(this) > a_Q
  */
  bool operator>(const int a_Q) const;

  /*!
    @brief Comparison operator. This uses lexicographical ordering
  */
  bool operator<(const MultiIndex& a_index) const;

  /*!
    @brief Comparison operator. Returns true if abs(this) < a_Q
  */
  bool operator<(const int a_Q) const;

  /*!
    @brief Comparison operator
  */
  bool operator <= (const int a_Q) const;

  /*!
    @brief Get next index in the permutation sequence
    @param[in] a_Q
  */
  void next(const int a_Q);

  /*!
    @brief Check if there are more lexigraphical permutations of the current index
  */
  bool hasNext(const int a_Q);

  /*!
    @brief Length of a multiindex 
  */
  static int abs(const MultiIndex& a_index);

  /*!
    @brief Power function for multi-indices
  */
  static Real pow(const RealVect& a_vec, const MultiIndex& a_index);
  
protected:

  /*!
    @brief Index
  */
  IntVect m_index;

  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Base factorial function
  */
  int factorial(const int a_n) const;
};

#endif
