#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MISC_H_
#define _MISC_H_

#include "CH_assert.H"
#include "REAL.H"
#include <cmath>
#include "BaseNamespaceHeader.H"

/**
   Returns a copy of the minimum of the two values.  Will work with
   any type T that supplies a T::operator<(const T, const T), or an
   equivalent construct, as well as a copy constructor.   In
   particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T> inline T Min(const T& a_a,
                                const T& a_b)
{
    return (a_a < a_b) ? a_a : a_b;
}

/**
   Returns a copy of the maximum of the two values.  Will work with
   any type T that supplies a T::operator>(const T, const T), or an
   equivalent construct, as well as a copy constructor.   In
   particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T> inline T Max (const T& a_a,
                                 const T& a_b)
{
    return (a_a > a_b) ? a_a : a_b;
}

/**
   Returns a copy of the absolute value of the value.  Will work with
   any type T that can be compared against zero with
   T::operator>(const T, const T), or an equivalent construct, as
   well as a copy constructor, and an T::operator-() or equivalent.
   In particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T> inline T Abs(const T& a_a)
{
    return (a_a > 0) ? a_a : -a_a;
}

/**
   Swaps the two values.  Type T must have a copy constructor and
   an assignment operator.
*/
template <class T> inline void Swap(T& a_a,
                                    T& a_b)
{
    T tmp = a_a;
    a_a   = a_b;
    a_b   = tmp;
}

//inline int pow(int a, int b)
//{
//  CH_assert(a>=0);
//  CH_assert(b>=0);
//  int rtn = 1;
//  for (; b>0; b--)
//    rtn*=a;
//  return rtn;
//}

inline int ipow(int a, int b)
{
  CH_assert(b>=0);
  int rtn = 1;
  for (; b>0; b--)
  {
     rtn *= a;
  }
  return rtn;
 // return std::pow(a,b);
}

inline Real ipow(const Real& a, const int& b)
{
  return std::pow(a, b);

//   double rtn = 1.0;
//   int r=Abs(b);
//   for (;r>0; r--)
//     rtn*=a;
//   return (b>0) ? rtn : 1.0/rtn;
}

/*------------------------------------------------------------------------------
 * The following rarely used functions are in namespace Misc to avoid conficts
 * and help identify their origin
 *----------------------------------------------------------------------------*/

namespace Misc
{

/// Test if integer is a power of 2
inline bool isPower2(const int a_i)
{
  if (a_i <= 0) return false;      // Catch <= 0
  unsigned i = a_i;
  while (!(i & 1))                 // Find first bit
    {
      i >>= 1;
    }
  return (i >> 1) ? false : true;  // Can't have any other bits
}

/// Test for class type
/**
   Use the SFINAE principle as described in 15.2.2 "Determining
   Class Types" in Vandevoorde and Josuttis "C++ Templates" book to
   see if T is class type.  E.g.,
     TypeTr<int>::IsClass == false
     TypeTr<Box>::IsClass == true
*/
template <typename T>
class TypeTr
{
  private:
  typedef char One;
  typedef struct
  {
    char a[2];
  } Two;
  template <typename C> static One test(int C::*);
  template <typename C> static Two test(...);
  public:
  enum
  {
    IsClass = sizeof(TypeTr<T>::template test<T>(0)) == 1
  };
};

}  // namespace Misc

#include "BaseNamespaceFooter.H"

#endif
