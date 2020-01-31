#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERMAPPEDCFINTERPSTENCIL_H_
#define _FOURTHORDERMAPPEDCFINTERPSTENCIL_H_

#include "FArrayBox.H"
#include "IntVect.H"
#include "Vector.H"

#include "NamespaceHeader.H"

/// Class to manage coarse-to-fine spatial interpolation to all fine cells within a single particular coarse cell

/**
 */
class FourthOrderMappedCFInterpStencil
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderMappedCFInterpStencil();

  /// Full constructor
  FourthOrderMappedCFInterpStencil(/// in range -(a_degree-1):(a_degree-1), specifies displacement to boundary in each dimension
                                   const IntVect&   a_bdryOffset,
                                   /// max degree of polynomials
                                   int              a_degree = 3,
                                   /// dimensions that are fixed, not interpolated
                                   Interval         a_fixedDims = Interval() )
  {
    define(a_bdryOffset, a_degree, a_fixedDims);
  }

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~FourthOrderMappedCFInterpStencil();

  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// in range -(a_degree-1):(a_degree-1), specifies displacement to boundary in each dimension
              const IntVect&   a_bdryOffset,
              /// max degree of polynomials
              int              a_degree = 3,
              /// dimensions that are fixed, not interpolated
              Interval         a_fixedDims = Interval() );

  const Vector<int>& indices() const
  { return m_indices; }

  const Vector<int>& indicesFull() const
  { return m_indicesFull; }

  const Vector<IntVect>& cells() const
  { return m_stencilCells; }

  const Vector<Real>& avgPowers() const
  { return m_avgPowers; }

protected:

  /// max degree of polynomials
  int m_degree;

  /// dimensions that are fixed, not interpolated
  Interval m_fixedDims;

  /// number of dimensions that are interpolated: SpaceDim - m_fixedDims.size()
  int m_nInterpDims;

  /// displacement from domain boundary in each dimension, in range -(a_degree-1):(a_degree-1)
  IntVect m_bdryOffset;

  /// number of coarse cells in the stencil
  int m_stencilSize;

  /// number of nonzero powers
  int m_numNonzeroPowers;

  /// length m_stencilSize, lists stencil cells
  Vector<IntVect> m_stencilCells;

  /// length m_stencilSize*(SpaceDim-m_fixedDims.size()), lists indices of stencil cells, grouped by cell
  Vector<int> m_indices;

  /// length m_stencilSize*SpaceDim, lists indices of stencil cells, grouped by cell
  Vector<int> m_indicesFull;

  /// length m_stencilSize*m_numNonzeroPowers, lists averages of <(x-x0)^p>j for stencil cells j, grouped by cell
  Vector<Real> m_avgPowers;

  /// whether define() has been called
  bool m_defined;

private:

  // Copy and assignment not allowed
  void operator=(const FourthOrderMappedCFInterpStencil& a_input);
  FourthOrderMappedCFInterpStencil(const FourthOrderMappedCFInterpStencil& a_input);
};

#include "NamespaceFooter.H"

#endif
