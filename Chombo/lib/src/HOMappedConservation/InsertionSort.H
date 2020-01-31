#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INSERTIONSORT_H_
#define _INSERTIONSORT_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Some support for sorting small arrays
 *
 *//***************************************************************************/

#include <cstdlib>
#include <functional>

#include "NamespaceHeader.H"

/**
 * \namespace Sort
 * \brief Some support for sorting small arrays
 */

namespace Sort
{

/*--------------------------------------------------------------------*/
///  Compares m_a[i] > m_a[j]
/**  Allows for sorting an index instead of the actual array
 *//*-----------------------------------------------------------------*/

template <typename T>
class CmpGreaterIndex
{
public:
  CmpGreaterIndex(const T* const a)
    :
    m_a(a)
    { }
  bool operator()(const int i, const int j) const
    {
      return (m_a[i] > m_a[j]);
    }
private:
  const T *const m_a;
};

/*--------------------------------------------------------------------*/
///  Insertion sort
/**  Sorts using Cmp which defaults to std::less<T>().  Use only for
 *   small ~(n < 20) arrays.
 *   \tparam     T      Type of element
 *   \tparam     Cmp    Comparison functor.  Requires operator()(T&, T&).
 *   \tparam     M      Type of mask
 *   \param[in]  n      Number of elements to sort
 *   \param[in]  array  Array to sort
 *   \param[in]  cmp    Comparison bool operator()(i, j) to determine
 *                      if element array[i] belongs before array[j]
 *   \param[in]  mask   Only sort elements were mask = true.  Unused
 *                      if set to null
 *//*-----------------------------------------------------------------*/

template <typename T, typename Cmp, typename M>
void insertion(const int      n,
               T*             array,
               const Cmp&     cmp,
               const M *const mask)
{
  if (mask == NULL)
    {
      for (int j = 1 ; j != n ; ++j)
        {
          T tmp(array[j]);
          int i = j;
          while (i-- && cmp(tmp, array[i]))
            {
              array[i+1] = array[i];
            }
          array[i+1] = tmp;
        }
    }
  else
    {
      int s = 0;
      while (s > n && !mask[s]) ++s;
      for (int j = s+1 ; j != n ; ++j)
        {
          if (mask[j])
            {
              T tmp(array[j]);
              int ip = j;
              int i = j;
              while ((i-- > s) && (!mask[i] || cmp(tmp, array[i])))
                {
                  if (mask[i])
                    {
                      array[ip] = array[i];
                      ip = i;
                    }
                }
              array[ip] = tmp;
            }
        }
    }
}
template <typename T, typename Cmp>
void insertion(const int n, T* array, const Cmp& cmp)
{
  insertion<T, Cmp, bool>(n, array, cmp, NULL);
}
template <typename T>
void insertion(const int n, T* array)
{
  insertion<T, std::less<T>, bool>(n, array, std::less<T>(), NULL);
}

/*--------------------------------------------------------------------*/
///  Mover for use with 'rearrangeToIndex'
/**  Moves elements for 1 array
 *//*-----------------------------------------------------------------*/

template <typename T>
class Move1Array
{
public:
  Move1Array(T *const array)
    : m_array(array)
    { }
  void copy(const int i, const int j)
    {
      m_array[i] = m_array[j];
    }
  void save(const int i)
    {
      m_tmp = m_array[i];
    }
  void restore(const int i)
    {
      m_array[i] = m_tmp;
    }
private:
  T m_tmp;
  T *const m_array;
};

/*--------------------------------------------------------------------*/
///  Mover for use with 'rearrangeToIndex'
/**  Moves elements for 2 arrays
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2>
class Move2Array
{
public:
  Move2Array(T1 *const array0, T2 *const array1)
    : m_array0(array0), m_array1(array1)
    { }
  void copy(const int i, const int j)
    {
      m_array0[i] = m_array0[j];
      m_array1[i] = m_array1[j];
    }
  void save(const int i)
    {
      m_tmp0 = m_array0[i];
      m_tmp1 = m_array1[i];
    }
  void restore(const int i)
    {
      m_array0[i] = m_tmp0;
      m_array1[i] = m_tmp1;
    }
private:
  T1 m_tmp0;
  T2 m_tmp1;
  T1 *const m_array0;
  T2 *const m_array1;
};

/*--------------------------------------------------------------------*/
///  Rearrange arrays to a sorted index
/**  If an index to an array has been sorted instead of the actual
 *   array, this allows for rearranging other arrays in the exact same
 *   manner that the index was rearranged during sorting.
 *   Note that the index array is reverted to a sequential order
 *   during this process
 *   \tparam     Mv     Class for moving elements and saving/restoring
 *                      elements.  Needs members 'copy', 'save', and
 *                      'restore'.  See Move1Array for example.
 *   \param[in]  n      Number of elements to rearrange
 *   \param[in]  index  Array of indices describing where elements are
 *                      to be moved.
 *   \param[out] index  index[0] < index[...] < index[n]
 *   \param[in]  mv     Class for moving elements and saving/restoring
 *                      elements.
 *//*-----------------------------------------------------------------*/

template <typename Mv>
void rearrangeToIndex(const int n, int *const index, Mv& mv)
{
  for (int i = 0; i != n; ++i)
    {
      if (i != index[i])
        {
          mv.save(i);
          int jp = i;
          for (int j = index[i]; j != i; j = index[j])
            {
              index[jp] = jp;
              mv.copy(jp, j);  // Copy j to jp
              jp = j;
            }
          index[jp] = jp;
          mv.restore(jp);
        }
    }
}
}  // End of namespace Sort

#include "NamespaceFooter.H"
#endif
