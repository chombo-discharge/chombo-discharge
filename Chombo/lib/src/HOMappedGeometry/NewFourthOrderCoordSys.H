#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWFOURTHORDERCOORDSYS_H_
#define _NEWFOURTHORDERCOORDSYS_H_

#include "REAL.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeQuadrature.H"
#include "NewCoordSys.H"

#include "NamespaceHeader.H"

/// Virtual base class encapsulating mapped-grid coordinate systems
/**
   The NewFourthOrderCoordSys class is a pure-virtual base class
   an API for performing operations in a mapped grid coordinate space.
   Implementations of specific coordinate systems will be provided by
   derived classes.

   The following pure virtual functions of NewCoordSys are defined in
   NewFourthOrderCoordSys.cpp :


  virtual void getNodeRealCoordinates(FArrayBox& a_nodeCoords,
                                      const Box& a_box) const;
  // calls virtual RealVect realCoord(const RealVect& a_Xi) const = 0;

  virtual int getNumN() const
  {
    return SpaceDim*SpaceDim;
  }

  virtual int getNcomponent(const int a_row, const int a_col) const
  {return (SpaceDim*a_col) + a_row;}

  virtual void cellVol(FArrayBox& a_vol,
                       const FluxBox& a_N,
                       const Box& a_box) const;
  // calls virtual RealVect realCoord(const RealVect& a_Xi) const = 0
  // on corners, then calls fourthOrderAverageFace and computeDivergence.

  virtual void getN(FluxBox& a_N, const Box& a_box) const;

  virtual void getAvgJ(FArrayBox& a_avgJ, const Box& a_box) const;
  // calls cellVol.

  virtual void getAvgJinverse(FluxBox& a_avgJinverse,
                              const Box& a_box) const;

  virtual void pointwiseJ(FArrayBox& a_J,
                          const FArrayBox& a_Xi,
                          const Box& a_box) const;
  // calls dXdXi.

  virtual void computeMetricTermProductAverage(FluxBox& a_product,
                                               const FluxBox& a_F,
                                               const FluxBox& a_N,
                                               const FluxBox& a_FforGrad,
                                               const Box& a_box,
                                               bool a_fourthOrder=true) const;

  virtual const RealVect& dx() const
  {
    return m_dx;
  }
*/
class NewFourthOrderCoordSys : public NewCoordSys
{
public:
  /// default constructor
  NewFourthOrderCoordSys();

  /**
     Destructor.
   */
  virtual ~NewFourthOrderCoordSys();

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

  /// return Cartesian XYZ locations of nodes
  /** nodeCoords should have dimension returned by dimension()
   */
  virtual void getNodeRealCoordinates(FArrayBox& a_nodeCoords,
                                      const Box& a_box) const;

  /// note that a_Xi is in mapped space.
  virtual Real  dXdXi(const RealVect& a_Xi,
                      int a_dirX,
                      int a_dirXi) const = 0;

  /// note that a_X is in mapped space.
  /** fills the destComp component of a_dxdXi
      with the derivative of x w/ respect to Xi
      in the dirX direction
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
  virtual int getNumN() const
  {
    return SpaceDim*SpaceDim;
  }

  /// index function into face-centered metric terms
  /** returns which component of the face-centered metric term in
      which N^row_col is stored.
      \note
      <ul>
        <li> \f$N\f$ is stored in column-major order (fortran).
        <li> Then using the same format, \f$N^T\f$ is stored in
             row-major order (C).
        <li> so (0,x) gives the start of a column of \f$N\f$ and
             a row of \f$N^T\f$.
        <li> If you want the computational flux in direction 'd',
             from the physical flux, 'f', take the scalar product
             of N(getNcomponent(0, d)) dot f
      </ul>
  */
  virtual int getNcomponent(const int a_row, const int a_col) const
  {return (SpaceDim*a_col) + a_row;}

  /// Same as getNcomponent but returns an IntVect for an 'a_col' (direction)
  static const IntVect& metricsCompDir(const int a_iDir)
    {
      static const IntVect metricsComp[SpaceDim] =
      {
        D_DECL6(IntVect(D_DECL6(0*SpaceDim + 0, 0*SpaceDim + 1,
                                0*SpaceDim + 2, 0*SpaceDim + 3,
                                0*SpaceDim + 4, 0*SpaceDim + 5)),
                IntVect(D_DECL6(1*SpaceDim + 0, 1*SpaceDim + 1,
                                1*SpaceDim + 2, 1*SpaceDim + 3,
                                1*SpaceDim + 4, 1*SpaceDim + 5)),
                IntVect(D_DECL6(2*SpaceDim + 0, 2*SpaceDim + 1,
                                2*SpaceDim + 2, 2*SpaceDim + 3,
                                2*SpaceDim + 4, 2*SpaceDim + 5)),
                IntVect(D_DECL6(3*SpaceDim + 0, 3*SpaceDim + 1,
                                3*SpaceDim + 2, 3*SpaceDim + 3,
                                3*SpaceDim + 4, 3*SpaceDim + 5)),
                IntVect(D_DECL6(4*SpaceDim + 0, 4*SpaceDim + 1,
                                4*SpaceDim + 2, 4*SpaceDim + 3,
                                4*SpaceDim + 4, 4*SpaceDim + 5)),
                IntVect(D_DECL6(5*SpaceDim + 0, 5*SpaceDim + 1,
                                5*SpaceDim + 2, 5*SpaceDim + 3,
                                5*SpaceDim + 4, 5*SpaceDim + 5)))
      };
      return metricsComp[a_iDir];
    }

  /// computes the volume flux on the faces
  virtual void volFlux(FluxBox& a_volFlux,
                       const FluxBox& a_Nt,
                       const Box& a_box) const;

  /// computes cell volumes
  /**
     Return cell volume in a_vol on cells of a_box,
     using a_N from getN() defined on a_box grown by 1.
   */
  virtual void cellVol(FArrayBox& a_vol,
                       const FluxBox& a_N,
                       const Box& a_box) const;

  /// computes integral of N over each face of a_box
  virtual void getN(FluxBox& a_N, const Box& a_box) const;

  /// computes cell-averaged J
  virtual void getAvgJ(FArrayBox& a_avgJ,
                       const FluxBox& a_volFlux,
                       const Box& a_box) const;

  /// computes cell-averaged J
  virtual void getAvgJ(FArrayBox& a_avgJ, const Box& a_box) const;

  /// computes cell-averaged 1/J
  virtual void getAvgJinverse(FluxBox& a_avgJinverse,
                              const Box& a_box) const;

  /// computes cell-averaged 1/J, using pre-computed cell-avg J
  virtual void getAvgJinverse(FluxBox& a_avgJinverse,
                              const FArrayBox& a_avgJ,
                              const Box& a_box) const;

  /// Jacobian evaluated at location Xi in mapped space
  virtual Real pointwiseJ(const RealVect& a_Xi) const;

  /// Jacobian evaluated at locations Xi in mapped space
  virtual void pointwiseJ(FArrayBox& a_J,
                          const FArrayBox& a_Xi,
                          const Box& a_box) const;

  /// Computes 4th-order average of product = N^T*F.  F can be a flux dyad
  virtual void computeMetricTermProductAverage(
    FluxBox&                   a_product,
    const FluxBox&             a_F,
    const FluxBox&             a_N,
    const int                  a_NNumComp,
    const FluxBox&             a_FforGrad,
    const Box&                 a_box,
    bool                       a_fourthOrder = true,
    Interval                   a_varIntervalProduct = Interval(),
    Interval                   a_varIntervalF = Interval(),
    int                        a_fluxSpaceContiguous = 1,
    const ProblemDomain *const a_problemDomainPtr = NULL) const;

  /// computes 4th-order average of product = N^T*F.  F is vector of a scalar
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
                                               bool a_fourthOrder=true) const;


  /// computes 4th-order average of product = N^T*F.  F is vector of a scalar
  /** if a_fourthOrder is false, then only do a second-order dot product;
      returns results in a_product on a_box;
      both a_F and a_N (from getN()) must be defined on a_box grown by 1.
   */
  virtual void computeMetricTermProductAverage(FluxBox& a_product,
                                               const FluxBox& a_F,
                                               const FluxBox& a_N,
                                               const Box& a_box,
                                               bool a_fourthOrder=true) const;

  /// Computes magnitude of <N^T> on each face
  //  This is not virtual because I don't think it should even be in this class
  //  -- see comments in NewFourthOrderCoordSys.cpp.
  void magnitudeN(FluxBox&       a_NMag,
                  const FluxBox& a_N,
                  const Box&     a_box) const;

  /// given N and Jinverse, computes N/J
  /** provided for convenience, since the multicomponent N/J product
      requires a bit of processing to work with the existing
      fourth-order multiplication function
  */
  virtual void computeNJinverse(FluxBox& a_NJinverse,
                                const FluxBox& a_Jinverse,
                                const FluxBox& a_N,
                                const Box& a_cellBox) const;


  /// access function to simplify things -- returns mapped-space cell spacing
  virtual const RealVect& dx() const
  {
    return m_dx;
  }

  /// this evaluates the script N values from equation 12 in Phil's notes
  /** note that a_Xi is in mapped space.
   */
  virtual Real getN(const RealVect& a_Xi, int a_s, int a_d, int a_d1) const;

  /// note that a_Xi is in mapped space.
  virtual Real getNMatrixEntry(const RealVect& a_Xi,
                               int a_s, int a_d, int a_d1,
                               int a_row, int a_column) const;

  ///  Integrates \f$\mathcal{N}^s\f$ on a specific set of codimension 2
  ///  hyperedges that are orthogonal to the given directions.
  virtual void integrateScriptN(FArrayBox& a_scrN,
                                const int a_dir0,
                                const int a_dir1,
                                const Box& a_box) const;

  virtual void incrementFaceMetricWithEdgeTerm(FArrayBox& a_faceMetrics,
                                               int a_faceDir,
                                               int a_edgeDir,
                                               const Box& a_box,
                                               const Side::LoHiSide& a_side) const;

  virtual void computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                           const Box& a_box,
                                           int a_faceDir,
                                           int a_dDir) const;


  virtual void computeTangentialGrad(FluxBox& a_gradPhi,
                                     const FluxBox& a_phiFace,
                                     const Box& a_box) const;


  virtual int tanGradComp(const int a_faceDir, const int a_tanDir,
                          const int a_comp) const;

  /// return row a_idir of contravariant metric
  virtual void contravariantMetric(FArrayBox& a_metric,
                                   int a_dir) const;

  /// convert vector components (in a_vectorIntv) of a_flux to orthonormal frame to that of this block
  virtual void orthonormalize(FluxBox& a_flux,
                              const Interval& a_vectorIntv) const;

  /// overwrite vector components (in a_vectorIntv) of a_fluxFab with orthonormalized
  virtual void orthonormalize(FArrayBox& a_fluxFab,
                              const Interval& a_vectorIntv) const;

  /// convert components in a_csComps of vector a_csFab in coordinate-system basis to components in a_orthoComps of vector a_orthoFab in orthonormal basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void orthonormalize(
                              const FArrayBox& a_csFab,
                              FArrayBox& a_orthoFab,
                              const Box& a_box,
                              int a_idir,
                              const IntVect& a_csComps,
                              const IntVect& a_orthoComps) const;

  /// convert vector components (in a_vectorIntv) of a_flux from orthonormal frame to that of this block
  virtual void deorthonormalize(FluxBox& a_flux,
                                const Interval& a_vectorIntv) const;

  /// overwrite vector components (in a_vectorIntv) of a_fluxFab with deorthonormalized
  virtual void deorthonormalize(FArrayBox& a_fluxFab,
                                const Interval& a_vectorIntv) const;

  /// convert components in a_orthoComps of vector a_orthoFab in orthonormal basis to components in a_csComps of vector a_csFab in coordinate-system basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void deorthonormalize(
                                const FArrayBox& a_orthoFab,
                                FArrayBox& a_csFab,
                                const Box& a_box,
                                int a_idir,
                                const IntVect& a_orthoComps,
                                const IntVect& a_csComps) const;
protected:

  RealVect m_dx;

  // which component directions to use when computing cell volumes
  Interval m_volInterval;

  EdgeQuadrature* m_quadraturePtr;
};

typedef NewFourthOrderCoordSys BlockCoordSys;


#include "NamespaceFooter.H"

#endif
