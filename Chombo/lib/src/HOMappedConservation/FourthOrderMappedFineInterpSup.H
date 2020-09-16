#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERMAPPEDFINEINTERPSUP_H_
#define _FOURTHORDERMAPPEDFINEINTERPSUP_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Supporting classes and routines for FourthOrderMappedFineInterp
 *
 *//***************************************************************************/

#include <cstring>

#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Routines
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
///  Calculate a binomial coefficient
/**  \param[in]         n
 *   \param[in]         k
 *   \return            Binomial coefficient
 *//*-----------------------------------------------------------------*/

inline int
binomial(const int n, int k)
{
  CH_assert(k >= 0);
  CH_assert(n >= k);
  if (k == 0) return 1;
  int cnum = n - k + 1;
  int cden = k;
  while (--k)
    {
      cnum *= n - k + 1;
      cden *= k;
    }
  return cnum/cden;
}

/*--------------------------------------------------------------------*/
///  Find the sequential index of a power
/**  The powers are indexed as follows (i.e for \f$d=3\f$ dimensions)  \verbatim
 *     int idx = 0
 *     for (px = 0; px <= m; ++px)
 *       for (py = 0; px+py <= m; ++py)
 *         for (pz = 0; px+py+pz <= m; ++py)
 *           ++idx                                                  \endverbatim
 *   where \f$m\f$ is the degree of the polynomial and we wish to find
 *   idx.  To compute the sequential index of any given power, we can
 *   use the relations
 *   \f[
 *      \mbox{number of powers} = {m + d \choose d}
 *   \f]
 *   and
 *   \f[
 *      \sum_{j=k}^n {j\choose k} = { n+1 \choose k+1 }
 *   \f]
 *   With these, the amount to add to the sequential index for a power
 *   at a spatial index is total number of powers remaining at this
 *   spatial index (remainder of \f$m\f$) minus the number of powers
 *   not used at this spatial index.  E.g, if \f$m=3\f$, \f$d=3\f$,
 *   and \f$p_x=1\f$, there are 2 powers left for the remaining 2
 *   dimensions,
 *   \f[
 *      {2+2 \choose 2}\,.
 *   \f]
 *   The increment to the sequential index is then
 *   \f[
 *      {2+3+1 \choose 2+1} - {2+2+1 \choose 2+1}
 *   \f]
 *   In general, this can be written for direction index \f$i\f$, in
 *   \f$d\f$ space dimensions, with \f$m_i\f$ giving the remaining
 *   available powers at \f$i\f$, and \f$p_i\f$ giving the power used
 *   at index \f$i\f$ as
 *   \f[
 *      {d-i+m_i \choose d-i} - {d-i+m_i-p_i \choose d-i } \,.
 *   \f]
 *
 *   \param[in]  a_m    Degree of the polynomial
 *   \param[in]  a_p    Power for each direction
 *   \return            Sequential index
 *//*-----------------------------------------------------------------*/

inline int
powerIndex(int a_m, IntVect a_p)
{
  int idx = 0;
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      const int remdir = SpaceDim - iDir;
      idx += binomial(remdir + a_m, remdir);
      a_m -= a_p[iDir];
      idx -= binomial(remdir + a_m, remdir);
    }
  return idx;
}


/*******************************************************************************
 */
///  Coordinate transformations
/**
 *   Provides for coordinate transformations.  The transformation is input in
 *   compact form where a vector of size SpaceDim describes for each direction,
 *   the corresponding direction (including sign) in transformed space.
 *   Directions are number from 1, e.g., [-2, -1, -3] transforms
 *   i -> -j, j -> -i, and k -> -k.  For transformation of IntVects, all
 *   coordinate systems are "fixed" at cell IntVect::Zero such that the
 *   transformation of IntVect::Zero is also IntVect::Zero.  See
 *   http://www.grc.nasa.gov/WWW/cgns/sids/cnct.html#GridConnectivity1to1
 *   for a more detailed explanation on transformations.
 *
 *   \note
 *   <ul>
 *     <li> The transformation description can be passed to Chombo Fortran as
 *          a 1-D array
 *          \verbatim CHF_CONST_I1D(ct.xctm(), 2*SpaceDim) \endverbatim
 *          and then used with a routine that duplicates
 *          CoordTransform::transform.
 *   </ul>
 *
 ******************************************************************************/

class CoordTransform
{
public:

  /// Default constructor
  CoordTransform();

  /// Constructor taking a compact description of the transformation
  CoordTransform(const int *const a_ctm);

  /// Copy constructor
  CoordTransform(const CoordTransform &a_ct);

  /// Assignment
  CoordTransform &operator=(const CoordTransform &a_ct);

  /// Defines the transformation by expanding a compact description
  void expand(const int *const a_ctm);

  /// Transforms an IntVect to the transformed space
  IntVect transform(const IntVect& a_iv) const;

  /// Returns the transformed direction associated with direction 'i'
  int indexAt(const int a_i) const;

  /// Returns the sign of the transformation associated with direction 'i'
  int signAt(const int a_i) const;

  /// Returns the expanded transformation description
  const int* xctm() const;

  /// Reverse the transformation to the other way
  void reverse();

private:
  /// Sign function
  int sign(const int a_i) const
    {
      return ( a_i < 0 ) ? -1 : 1;
    }

  int m_xctm[2*SpaceDim];             ///< An expanded description of the
                                      ///< transformation.  For direction i:
                                      ///< m_xctm[2*i]   = new direction
                                      ///< m_xctm[2*i+1] = sign
};


/*******************************************************************************
 *
 * Class CoordTransform: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//   Default sets m_xctm everywhere to zero.
/**  \param[in]  a_ctm  Compact description of the transformation
 *//*-----------------------------------------------------------------*/

inline
CoordTransform::CoordTransform()
{
  std::memset(m_xctm, 0, 2*SpaceDim*sizeof(int));
}

/*--------------------------------------------------------------------*/
//   Constructor taking a compact description of the transformation
/**  \param[in]  a_ctm  Compact description of the transformation
 *//*-----------------------------------------------------------------*/

inline
CoordTransform::CoordTransform(const int *const a_ctm)
{
  expand(a_ctm);
}

/*--------------------------------------------------------------------*/
//   Copy constructor
/**  \param[in]  a_ct   CoordTransform to copy
 *//*-----------------------------------------------------------------*/

inline
CoordTransform::CoordTransform(const CoordTransform &a_ct)
{
  std::memcpy(m_xctm, a_ct.m_xctm, 2*SpaceDim*sizeof(int));
}

/*--------------------------------------------------------------------*/
//   Assignment
/**  \param[in]  a_ct   CoordTransform to copy
 *   \return            Updated with assignment
 *//*-----------------------------------------------------------------*/

inline CoordTransform &
CoordTransform::operator=(const CoordTransform &a_ct)
{
  if (&a_ct != this)
    {
      std::memcpy(m_xctm, a_ct.m_xctm, 2*SpaceDim*sizeof(int));
    }
  return *this;
}

/*--------------------------------------------------------------------*/
//   Defines the transformation by expanding a compact description
/**  \param[in]  a_ctm  Compact description of the transformation
 *//*-----------------------------------------------------------------*/

inline void
CoordTransform::expand(const int *const a_ctm)
{
  for (int i = 0; i != SpaceDim; ++i)
    {
      m_xctm[2*i]   = std::abs(a_ctm[i]) - 1;
      m_xctm[2*i+1] = sign(a_ctm[i]);
    }
}

/*--------------------------------------------------------------------*/
//   Transforms an IntVect to the transformed space
/**  \param[in]  a_iv   IntVect to transform
 *   \return            Tranformed IntVect
 *//*-----------------------------------------------------------------*/

inline IntVect
CoordTransform::transform(const IntVect& a_iv) const
{
  IntVect tiv;
  D_EXPR6(tiv[m_xctm[ 0]] = m_xctm[ 1]*a_iv[0],
          tiv[m_xctm[ 2]] = m_xctm[ 3]*a_iv[1],
          tiv[m_xctm[ 4]] = m_xctm[ 5]*a_iv[2],
          tiv[m_xctm[ 6]] = m_xctm[ 7]*a_iv[3],
          tiv[m_xctm[ 8]] = m_xctm[ 9]*a_iv[4],
          tiv[m_xctm[10]] = m_xctm[11]*a_iv[5]);
  return tiv;
}

/*--------------------------------------------------------------------*/
//   Returns the transformed direction associated with direction 'i'
/**  \param[in]  a_i    Input direction
 *   \return            Transformed direction
 *//*-----------------------------------------------------------------*/

inline int
CoordTransform::indexAt(const int a_i) const
{
  return m_xctm[2*a_i];
}

/*--------------------------------------------------------------------*/
//   Returns the sign of the transformation associated with direction
//   'i'
/**  \param[in]  a_i    Input direction
 *   \return            Sign associated with transformation
 *//*-----------------------------------------------------------------*/

inline int
CoordTransform::signAt(const int a_i) const
{
  return m_xctm[2*a_i + 1];
}

/*--------------------------------------------------------------------*/
//   Returns the expanded transformation description
/**  \return            Exanded transformation decsription
 *//*-----------------------------------------------------------------*/

inline const int*
CoordTransform::xctm() const
{
  return m_xctm;
}

/*--------------------------------------------------------------------*/
//   Reverse the transformation to the other way
/*--------------------------------------------------------------------*/

inline void
CoordTransform::reverse()
{
  int origxctm[2*SpaceDim];
  std::memcpy(origxctm, m_xctm, 2*SpaceDim*sizeof(int));
  for (int i = 0; i != SpaceDim; ++i)
    {
      const int j = origxctm[2*i];
      m_xctm[2*j] = i;
      m_xctm[2*j+1] = origxctm[2*i+1];
    }
}

#include "NamespaceFooter.H"
#endif
