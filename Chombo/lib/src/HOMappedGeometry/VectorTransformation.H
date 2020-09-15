#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VECTORTRANSFORMATION_H_
#define _VECTORTRANSFORMATION_H_

#include "IntVect.H"
#include "RealVect.H"
#include "Tuple.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe transformation of SpaceDim-vectors, typically from one block to another.
class VectorTransformation
{

  friend class VectorAvgTransformation;

public:

  /// null constructor leaves object in undefined state
  VectorTransformation()
  { }

  /// set one component of transformation matrix
  void setComponent(int a_dirDst,
                    int a_dirSrc,
                    Real a_val)
  { m_coeffs[matrixIndex(a_dirDst, a_dirSrc)] = a_val; }

  Real component(int a_dirDst,
                 int a_dirSrc) const
  { return m_coeffs[matrixIndex(a_dirDst, a_dirSrc)]; }

  /// only for sending stuff to Fortran
  const Real* dataPtr() const
  { return m_coeffs; }

  /// apply this tranformation to vector
  RealVect transform(const RealVect& a_vec) const;

  /// increment a vector-transformed cell-centered FAB
  void addTransformFab(FArrayBox& a_dstFab,
                       const FArrayBox& a_srcFab,
                       const Box& a_bx,
                       Real a_weight = 1.) const;

  /// increment a vector-transformed cell-averaged FAB
  void addTransformFabAvg(FArrayBox& a_dstFab,
                          const FArrayBox& a_srcFab,
                          const FArrayBox& a_gradFab,
                          const Vector<VectorTransformation>& a_vecTransformGrad,
                          const Box& a_bx,
                          Real a_weight = 1.) const;
  ///
  /**
     Return the inverse transformation.
   */
  VectorTransformation inverse() const;

  ///
  /**
     Return the composite transformation:
     that is, this transformation followed by a_next.
   */
  VectorTransformation compose(const VectorTransformation& a_next) const;

  VectorTransformation& operator+= (const VectorTransformation& p);

  friend VectorTransformation operator+ (const VectorTransformation& s,
                                         const VectorTransformation& p);

  VectorTransformation& operator-= (const VectorTransformation& p);

  friend VectorTransformation operator- (const VectorTransformation& s,
                                         const VectorTransformation& p);

  VectorTransformation& operator*= (Real s);

  VectorTransformation operator* (Real s) const;

  friend VectorTransformation operator* (Real            s,
                                         const VectorTransformation& p);

  VectorTransformation& operator/= (Real s);

  VectorTransformation operator/ (Real s) const;

  friend VectorTransformation operator/ (Real            s,
                                         const VectorTransformation& p);
                                         

  ///
  /** The identity transformation, which leaves a vector the same.
   */
  static const VectorTransformation Identity;

  ///
  /** Initializes VectorTransformation::Identity .
   */
  static int InitStatics();

protected:

  const static int s_ncoeffs = SpaceDim*SpaceDim;

  // Tuple<RealVect, SpaceDim> m_matrix;
  // Here RealVect m_matrix[idirDst] is for one destination.
  // Length s_ncoeffs.
  Real m_coeffs[s_ncoeffs];

  // order in m_coeffs:  [idirDst,idirSrc] ==
  // [0,0], [0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2];
  // |-- m_matrix[0] --|  |-- m_matrix[1] --| |--- m_matrix[2] --|
  // hence index in m_coeffs is SpaceDim*idirDst + idirSrc.

  static inline int matrixIndex(int a_dirDst, int a_dirSrc)
  {
    return (SpaceDim*a_dirDst + a_dirSrc);
  }

  Vector< Vector<Real> > vecvec(const Real* a_coeffs) const;

  Real determinant(const Real* a_matrix) const;

  Real determinant(const Vector< Vector<Real> >& a_matrix) const;

  Real matrixMinor(const Vector< Vector<Real> >& a_matrix,
                   int a_row,
                   int a_col) const;
};

//
//// Static initialization.  Gotta make sure there are no static object
//// definitions above here (except possibly stuff in headers).  Otherwise,
//// the danger is that some static object's constructor calls
//// VectorTransformation::Identity, the very things the following
//// definition is supposed to initialize.
//// (I got this from IntVect.H)
static int s_dummyForVectorTransformationH = VectorTransformation::InitStatics();

#include "NamespaceFooter.H"

#endif // include guard
