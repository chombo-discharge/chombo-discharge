#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWCOORDSYS_H_
#define _NEWCOORDSYS_H_

#include "REAL.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "ProblemDomain.H"
#include "VectorTransformation.H"

#include "NamespaceHeader.H"

/// Virtual base class encapsulating mapped-grid coordinate systems
/**
   The CoordSys class is a pure-virtual base class providing an API
   for performing operations in a mapped grid coordinate space.
   Implementations of specific coordinate systems will be provided by
   derived classes.
*/
class NewCoordSys
{
public:
  /// default constructor
  NewCoordSys()
    : m_time(0.)
  {
  }

  /**
     Destructor.
   */
  virtual ~NewCoordSys()
  {;}

  /// given coordinate in mapped space, return its location in real space
  virtual RealVect realCoord(const RealVect& a_Xi) const = 0;

  /// given coordinate in real space, return its location in the mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const = 0;

  /// given coordinates in mapped space, return locations in real space
  virtual void realCoord(FArrayBox& a_x, const FArrayBox& a_Xi,
                         const Box& a_box) const;

  /// given coordinate in real space, return its location in the mapped space
  virtual void mappedCoord(FArrayBox& a_Xi, const FArrayBox& a_x,
                           const Box& a_box) const;

  /// given indices in a_box, return centers of cells (or of faces, or of edges) in mapped coordinates.  The function does the right thing based on the centering of a_box, which must be the same centering as a_Xi.
  virtual void getCenterMappedCoordinates(FArrayBox& a_Xi,
                                          const Box& a_box) const;

  /// given indices, return centers of cells in mapped coordinates
  virtual RealVect centerMappedCoordinates(const IntVect& a_iv) const;

  /// return Cartesian XYZ locations of nodes
  /** nodeCoords should have dimension returned by dimension()
   */
  virtual void getNodeRealCoordinates(FArrayBox& a_nodeCoords,
                                      const Box& a_box) const =0;


  /// note that a_Xi is in mapped space.
  virtual Real  dXdXi(const RealVect& a_Xi,
                      int a_dirX,
                      int a_dirXi) const = 0;

  /// note that a_Xi is in mapped space.
  /** fills the destComp component of a_dxdXi
      with the derivative of x w/ respect to Xi
      in the dirX direction
      Note that the default implementation of this in
      FourthOrderCoordSys class throws an error, since there is no
      way to get real=space coordinate x from index location except in
      the derived class..
      (default implementation is there in case derived class doesn't
      provide it)
  */
  virtual void dXdXi(FArrayBox& a_dxdXi,
                     const FArrayBox& a_Xi,
                     int a_destComp,
                     int a_dirX,
                     int a_dirXi,
                     const Box& a_box) const;


  /// returns number of components in the metric term matrix (N)
  virtual int getNumN() const = 0;

  /// index function into face-centered metric terms
  /** returns which component of the face-centered metric term in
      which N^row_col is stored
  */
  virtual int getNcomponent(const int a_row, const int a_col) const = 0;


  /// returns dimensionality of the space
  /**
      The dimensionality of the space may not be SpaceDim if, for example
      we're computing on a manifold. One example of this is the
      surface-of-a-sphere mapping, which is a 2D surface in a 3D space.
      For that case, dimension() would return 3 even though SpaceDim is 2.
      Default implementation returns SpaceDim.
  */
  virtual int dimension() const
  {
    return SpaceDim;
  }

  /// computes cell volumes
  /**
     Return cell volume in a_vol on cells of a_box,
     using a_N from getN() defined on a_box grown by 1.
   */
  virtual void cellVol(FArrayBox& a_vol,
                       const FluxBox& a_N,
                       const Box& a_box) const = 0;

  /// computes integral of N over each face of a_box
  virtual void getN(FluxBox& a_N, const Box& a_box) const = 0;

  /// computes cell-averaged J
  virtual void getAvgJ(FArrayBox& a_avgJ, const Box& a_box) const = 0;


  /// computes face-averaged 1/J
  virtual void getAvgJinverse(FluxBox& a_avgJinverse,
                              const Box& a_box) const = 0;
  /// Jacobian evaluated at location Xi in mapped space
  virtual Real pointwiseJ(const RealVect& a_Xi) const = 0;

  /// Jacobian evaluated at locations Xi in mapped space
  virtual void pointwiseJ(FArrayBox& a_J,
                          const FArrayBox& a_Xi,
                          const Box& a_box) const;

  /// returns integral of divergence over mapped-grid cells
  virtual void computeDivergence(FArrayBox& a_divF,
                                 const FluxBox& a_F,
                                 const FluxBox& a_N,
                                 const Box& a_box,
                                 Interval& divInterval) const;

  /// computes 4th-order average of product = N^T*F
  /** if a_fourthOrder is false, then only do a second-order dot product
      aFforGrad is the F used to compute grad_perp(F);
      returns results in a_product on a_box;
      both a_FforGrad and a_N (from getN()) must be defined on a_box grown by 1.
   */
  virtual void computeMetricTermProductAverage(FluxBox& a_product,
                                               const FluxBox& a_F,
                                               const FluxBox& a_N,
                                               const FluxBox& a_FforGrad,
                                               const Box& a_box,
                                               bool a_fourthOrder=true) const =0;


  /// computes 4th-order average of product = N^T*F
  /** if a_fourthOrder is false, then only do a second-order dot product;
      returns results in a_product on a_box;
      both a_F and a_N (from getN()) must be defined on a_box grown by 1.
   */
  virtual void computeMetricTermProductAverage(FluxBox& a_product,
                                               const FluxBox& a_F,
                                               const FluxBox& a_N,
                                               const Box& a_box,
                                               bool a_fourthOrder=true) const;

  /// return matrix that transforms vector components from mapped space to real space at point a_Xi in mapped space
  virtual VectorTransformation vectorTransformationMappedToReal(const RealVect& a_Xi) const;

  /// return matrix that transforms vector components from real space to mapped space at point a_Xi in mapped space
  virtual VectorTransformation vectorTransformationRealToMapped(const RealVect& a_Xi) const;

  /// return matrix that transforms vector components from mapped space to real space at center of cell a_iv
  virtual VectorTransformation vectorTransformationMappedToRealCenter(const IntVect& a_iv) const;

  /// return matrix that transforms vector components from real space to mapped space at center of cell a_iv
  virtual VectorTransformation vectorTransformationRealToMappedCenter(const IntVect& a_iv) const;

  /// transform a FAB of SpaceDim-vectors from mapped-coordinate basis to real-coordinate basis at cell centers
  virtual void vectorTransformMappedToRealCenterFab(FArrayBox& a_vectorFab) const;

  /// transform a FAB of SpaceDim-vectors from mapped-coordinate basis to real-coordinate basis, averaged over cells
  virtual void vectorTransformMappedToRealAvgFab(FArrayBox& a_vectorFab,
                                                 const FArrayBox& a_gradVectorFab) const;

  /// transform a FAB of SpaceDim-vectors from real-coordinate basis to mapped-coordinate basis at cell centers
  virtual void vectorTransformRealToMappedCenterFab(FArrayBox& a_vectorFab) const;

  /// transform a FAB of SpaceDim-vectors from real-coordinate basis to mapped-coordinate basis, averaged over cells
  virtual void vectorTransformRealToMappedAvgFab(FArrayBox& a_vectorFab,
                                                 const FArrayBox& a_gradVectorFab) const;

  /// access function to simplify things -- returns mapped-space cell spacing
  virtual const RealVect& dx() const =0;

  /// Set time for moving grids
  virtual void setTime(const Real& a_time) const;

protected:

  /// internal function to transform a FAB of SpaceDim-vectors from one basis to another, averaged over cells, given BaseFab of vector transformations at cell centers
  virtual void vectorTransformInternalAvgFab(FArrayBox& a_vectorFab,
                                             const FArrayBox& a_gradVectorFab,
                                             const BaseFab<VectorTransformation>& a_tfmCenterFab) const;

  mutable Real m_time;                ///< Const is already well entrenched
};


/*******************************************************************************
 *
 * Class NewCoordSys: inline member definitions
 *
 ******************************************************************************/

// Set time for moving grids
inline void
NewCoordSys::setTime(const Real& a_time) const
{
  m_time = a_time;
}

/// pure-virtual base class defining factory interface for CoordSys
class
NewCoordSysFactory
{


public:

  NewCoordSysFactory()
  {
  }

  virtual ~NewCoordSysFactory()
  {
  }

  virtual NewCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                   const RealVect& a_dx) const = 0;
};

#include "NamespaceFooter.H"

#endif
