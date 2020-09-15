#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SHAPEARRAY_H_
#define _SHAPEARRAY_H_

#include <cstddef>
#include <tuple>
#include "IntVect.H"

#include "NamespaceHeader.H"

#if CXXSTD>=14
namespace shape
{

/// Returns length of dimension D from a tuple of lengths
/** If highest dimension, returns 0
 *  \tparam IIx        Indexing type
 *  \tparam D          Dimension to query length
 *  \tparam Tuple      Tuple type
 *  \param[in]  a_lens Tuple of all dimension lengths
 *  \return            D = rank : 0
 *                     D < rank : length of the dimension
 */
template <typename IIx, std::size_t D, typename Tuple>
constexpr IIx getLen(const Tuple& a_lens) noexcept
{
  static_assert(std::is_integral<
                typename std::tuple_element<D, Tuple>::type
                >::value,
                "ShapeArray requires an integer dimension");
  return (D == std::tuple_size<Tuple>::value - 1) ? 0 : std::get<D>(a_lens);
}


/*******************************************************************************
 */
/// Shapes linear memory according to some rank and dimensions and gives access
/**
 *  Template recursion is used to store the size of each dimension and at D==1,
 *  a pointer to the data is also stored.  Access is provided through
 *  multidimensional notation x[i2][i1][i0] for rank = 3.  C (row) ordering is
 *  assumed, where i0 (D==1) is unit stride.  The size of the highest dimension
 *  (D==rank) is not required and this space is reused during access.
 *
 *  3-D (rank==3) example:
 *  We want an array of size i=4 x j=3 x k=2 where i is unit stride.
 *  Construct with a pointer to linear memory and the size of each dimension.
 *                                                                   \code{.cpp}
 *    double *data = new double[4*3*2];
 *    ShapeArray<double, 3> array3d(data, 4, 3, 2);
 *                                                                      \endcode
 *  Alternatively, use the helper function where template parameters are
 *  deduced.
 *                                                                   \code{.cpp}
 *    auto array3d = make_MD_array(data, 4, 3, 2);
 *                                                                      \endcode
 *  Access the array using operator[] for each dimension
 *                                                                   \code{.cpp}
 *    for (int k = 0; k != 2; ++k)
 *    for (int j = 0; j != 3; ++j)
 *    for (int i = 0; i != 4; ++i)
 *      {
 *        assert(&array3d[k][j][i] == &data[(k*3 + j)*4 + i]);
 *      }
 *                                                                      \endcode
 *
 *  Note that operator[] returns a ShapeArray object of rank D-1.  When D=1,
 *  the data is returned.  During access, the size of the highest dimension
 *  is modified to tabulate the math on the right-hand-side of the equality
 *  above; lower dimensions are propagated.
 *
 *  Performance:
 *  Compilers are very effective at optimizing the access.  Constructing the
 *  ShapeArray is costly, so you would not want to do that for every access.
 *  However access is equivalent in performance to built-in Fortran arrays or
 *  in C, casting to pointers-to-VLAs to make use of built-in C-style
 *  multidimensional indexing.
 *
 *  Non-zero lower bounds:
 *  The indexing assumes a zero-based lower bound.  If otherwise, you must
 *  offset the array data.  For a Chombo box, this can be done with.
 *                                                                   \code{.cpp}
 *    double *offsetData = data + box.index(IntVect_zero);
 *                                                                      \endcode
 *  If you ever shift the box, make sure you reset the data with a new offset.
 *
 *  \tparam T           Type of data in the array
 *  \tparam D           Rank of the array (during recursion, this is the
 *                      current dimension)
 *  \tparam IIx         Type of integer for indexing the array (default
 *                      std::ptrdiff_t)
 *  
 *//*+*************************************************************************/

template <typename T, std::size_t D, typename IIx = std::ptrdiff_t>
struct ShapeArray : ShapeArray<T, D-1, IIx>
{
  /// Default constructor
  ShapeArray() noexcept
    :
    ShapeArray<T, D-1, IIx>()
    { }

  /// Construct from integer sequence
  /** \tparam Is...       Integer sequence of array dimensions
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_lens  Size of each dimension
   */
  template <typename... Is>
  ShapeArray(T* a_data, Is... a_lens) noexcept
    :
    ShapeArray<T, D-1, IIx>(a_data, a_lens...),
    m_thisDlen(getLen<IIx, D-1>(std::make_tuple(a_lens...)))
    { }

  /// Construct from a static vector.  The rank of the array must be SpaceDim+1
  /** \tparam T2          Type of data in the vector
   *  \tparam D2          Size of the vector (must be >= SpaceDim)
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_lens  Size of each dimension
   */
  // template <typename T2, stc::array_size_type D2>
  // ShapeArray(T* a_data, const stc::Vector<T2, D2>& a_lens) noexcept
  //   :
  //   ShapeArray<T, D-1, IIx>(a_data, a_lens),
  //   m_thisDlen((D-1 == SpaceDim) ? 0 : a_lens[D-1])
  //   { }
  ShapeArray(T* a_data, const IntVect& a_lens) noexcept
    :
    ShapeArray<T, D-1, IIx>(a_data, a_lens),
    m_thisDlen((D-1 == SpaceDim) ? 0 : a_lens[D-1])
    { }

  /// Copy constructor
  ShapeArray(const ShapeArray& a_other) noexcept
    :
    ShapeArray<T, D-1, IIx>((ShapeArray<T, D-1, IIx>)a_other),
    m_thisDlen(a_other.m_thisDlen)
    { }

  /// Assignment operator
  ShapeArray& operator=(const ShapeArray& a_other) = default;

  /// Define from integer sequence
  /** \tparam Is...       Integer sequence of array dimensions
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_lens  Size of each dimension
   */
  template <typename... Is>
  void define(T* a_data, Is... a_lens) noexcept
    {
      ShapeArray<T, D-1, IIx>::define(a_data, a_lens...);
      m_thisDlen = getLen<IIx, D-1>(std::make_tuple(a_lens...));
    }

  /// Define from a static vector.  The rank of the array must be SpaceDim+1
  /** \tparam T2          Type of data in the vector
   *  \tparam D2          Size of the vector (must be >= SpaceDim)
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_lens  Size of each dimension
   */
  // template <typename T2, stc::array_size_type D2>
  // void define(T* a_data, const stc::Vector<T2, D2>& a_lens) noexcept
  //   {
  //     ShapeArray<T, D-1, IIx>::define(a_data, a_lens);
  //     m_thisDlen = (D-1 == SpaceDim) ? 0 : a_lens[D-1];
  //   }
  void define(T* a_data, const IntVect& a_lens) noexcept
    {
      ShapeArray<T, D-1, IIx>::define(a_data, a_lens);
      m_thisDlen = (D-1 == SpaceDim) ? 0 : a_lens[D-1];
    }

  /// Reset only the data
  /** This is often used if the bounds change (e.g., a box is shifted)
   *  \param[in]  a_data  Linear memory to reshaped
   */
  void resetData(T* a_data) noexcept
    {
      ShapeArray<T, 1, IIx>::resetData(a_data);
    }

  /// Make it obvious that the array is not defined
  void clear() noexcept
    {
      ShapeArray<T, D-1, IIx>::clear();
      m_thisDlen = 0;
    }

  /// Access the data
  /** \param[in]  a_idx   Index for this dimension
   *  \return             A temporary ShapeArray of one lower dimension
   */
  auto operator[](const IIx a_idx) const noexcept
    {
      return ShapeArray<T, D-1, IIx>(m_thisDlen + a_idx,
                                     // propagate the original array
                                     (ShapeArray<T, D-1, IIx>)(*this));
    }

  // Higher dimension object can access our data
  friend struct ShapeArray<T, D+1, IIx>;

protected:

  /// Access builder
  /** Higher-dimension indexing is stored in m_thisDlen.  Lower dimensions are
   *  propagated.
   */
  ShapeArray(const IIx a_idx, const ShapeArray& a_other) noexcept
    :
    ShapeArray<T, D-1, IIx>((ShapeArray<T, D-1, IIx>)a_other),
    m_thisDlen(a_idx*a_other.m_thisDlen)
    { }

  IIx m_thisDlen;                     ///< Length of this dimension (or storage
                                      ///< of index math if this is a temporary
                                      ///< during access)
}; 

// Specialization for D==1
template <typename T, typename IIx>
struct ShapeArray<T, 1, IIx>
{
  ShapeArray() noexcept
    :
    m_data(nullptr)
    { }
  template <typename... Is>
  ShapeArray(T* a_data, Is... a_lens) noexcept
    :
    m_thisDlen(getLen<IIx, 0>(std::make_tuple(a_lens...))),
    m_data(a_data)
    { }
  // template <typename T2, stc::array_size_type D2>
  // ShapeArray(T* a_data, const stc::Vector<T2, D2>& a_lens) noexcept
  //   :
  //   m_thisDlen(a_lens[0]),
  //   m_data(a_data)
  //   {
  //     static_assert(D2 >= SpaceDim, "ShapeArray constructed from static vector "
  //                   "with insufficient size");
  //   }
  ShapeArray(T* a_data, const IntVect& a_lens) noexcept
    :
    m_thisDlen(a_lens[0]),
    m_data(a_data)
    { }
  ShapeArray(const ShapeArray& a_other) noexcept
    :
    m_thisDlen(a_other.m_thisDlen),
    m_data(a_other.m_data)
    { }
  ShapeArray& operator=(const ShapeArray& a_other) = default;
  template <typename... Is>
  void define(T* a_data, Is... a_lens) noexcept
    {
      m_thisDlen = getLen<IIx, 0>(std::make_tuple(a_lens...));
      m_data = a_data;
    }
  // template <typename T2, stc::array_size_type D2>
  // void define(T* a_data, const stc::Vector<T2, D2>& a_lens) noexcept
  //   {
  //     static_assert(D2 >= SpaceDim, "ShapeArray defined from static vector "
  //                   "with insufficient size");
  //     m_thisDlen = a_lens[0];
  //     m_data = a_data;
  //   }
  void define(T* a_data, const IntVect& a_lens) noexcept
    {
      m_thisDlen = a_lens[0];
      m_data = a_data;
    }
  void resetData(T* a_data) noexcept
    {
      m_data = a_data;
    }
  void clear() noexcept
    {
      m_thisDlen = 0;
      m_data = nullptr;
    }

  /// Access the data
  /** \param[in]  a_idx   Index for this dimension
   *  \return             Modifiable data at requested element
   */
  T& operator[](const IIx a_idx) const noexcept
    {
      return *(m_data + m_thisDlen + a_idx);
    }

  /// Get pointer to the 1-D data
  T* data() const noexcept
    {
      return m_data + m_thisDlen;
    }

  // Higher dimension object can access our data
  friend struct ShapeArray<T, 2, IIx>;

protected:

  /// Access builder
  /** Higher-dimension indexing is stored in m_thisDlen.
   */
  ShapeArray(const IIx a_idx, const ShapeArray& a_other) noexcept
    :
    m_thisDlen(a_idx*a_other.m_thisDlen),
    m_data(a_other.m_data)
    { }

  IIx m_thisDlen;                     ///< Length of this dimension
  T* m_data;                          ///< Linear data
};

/// Helper function that builds a ShapeArray deducing all template parameters
/** \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_lens  Size of each dimension
 *  \return     ShapeArray with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Is>
inline auto make_MD_array(T* a_data, Is... a_lens) noexcept
{
  return ShapeArray<T, sizeof...(a_lens), int>(a_data, a_lens...);
}

}  // namespace shape
#endif

#include "NamespaceFooter.H"

#endif  /* ! defined _SHAPEARRAY_H_ */
