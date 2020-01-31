#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERE2DPANELCS_H_
#define _CUBEDSPHERE2DPANELCS_H_

#include <utility>
using std::pair;

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "RealVect.H"
#include "IntVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "Vector.H"
#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

/// Basic single-panel cubed sphere Coordinates
/**
   The CubedSphere2DPanelCS class implements the NewFourthOrderCoordSys class
   for cubed sphere panels in equiangular coordinates.
*/
class CubedSphere2DPanelCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     nPanel
       (0-3 Equatorial, 4 North Pole, 5 South Pole) and defines the
       panel on the cubed sphere described by this object.
     dDX
       The cell spacing on this panel.
     ix
       Index of the bottom-left edge of this panel in index space.
  */
  CubedSphere2DPanelCS(int nPanel, RealVect& dDX, IntVect& ix);

  void setFlatMap(bool a_flatMap);

  /**
     Destructor.
   */
  virtual ~CubedSphere2DPanelCS();

  //void panel(int nPanel) {m_nPanel = nPanel;}

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  void realCoord(FArrayBox& a_x,
                 const FArrayBox& a_Xi,
                 const Box& a_box) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  /// return Cartesian XYZ locations of cell centers
  virtual void getNodeRealCoordinates(FArrayBox& a_nodeCoords,
                                      const Box& a_box) const;

  // return cell center positions in mapped coordinates
  virtual void getCellMappedCoordinates(FArrayBox& a_cellCoords,
                                        const Box& a_box) const;

  // return face center positions in mapped coordinates
  virtual void getFaceMappedCoordinates(FArrayBox& a_faceCoords,
                                        const int a_dir,
                                        const Box& a_box) const;

  // return node positions in mapped coordinates
  virtual void getNodeMappedCoordinates(FArrayBox& a_nodeCoords,
                                        const Box& a_box) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  /// returns dimensionality of the space (3, or 2 if flat)
  /**
  */
  virtual int dimension() const
  {
    return m_realDim;
  }

  /// computes cell volumes
  virtual void cellVol(FArrayBox& a_vol,
                       const FluxBox& a_N,
                       const Box& a_box) const;

  /// computes NT * F.
  virtual void volFlux(FluxBox& a_volFlux,
                       const FluxBox& a_Nt,
                       const Box& a_box) const;

  /// computes face-centered N
  virtual void getN(FluxBox& a_N, const Box& a_box) const;

  /// computes cell-averaged J
  virtual void getAvgJ(FArrayBox& a_avgJ, const Box& a_box) const;

  /// computes cell-averaged J; ignores a_volFlux
  virtual void getAvgJ(FArrayBox& a_avgJ,
                       const FluxBox& a_volFlux,
                       const Box& a_box) const;

  /// computes face-averaged 1/J
  virtual void getAvgJinverse(FluxBox& a_avgJinverse,
                              const Box& a_box) const;

  /// Jacobian evaluated at location X in real space
  virtual Real pointwiseJ(const RealVect& a_X) const;

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

  /// this evaluates the script N values from equation 12 in Phil's notes
  /** note that a_Xi is in mapped space.
   */
  virtual Real getN(const RealVect& a_Xi, int a_s, int a_d, int a_d1) const;

  virtual Real getNMatrixEntry(const RealVect& a_Xi,
                               int a_s, int a_d, int a_d1,
                               int a_row, int a_column) const;

  /// transform a point from mapped-coordinate basis to longitude-latitude
  /// coordinate basis
  void pointTransformEquiangularToLonLat(const RealVect& a_xi,
                                         RealVect & a_rllXi) const;

  /// transform a point from mapped-coordinate basis to Cartesian
  /// coordinate basis
  void pointTransformEquiangularToCartesian(const RealVect& a_xi,
                                            Real * a_xyz) const;

  /// transform a FAB of points from mapped-coordinate basis to latitude-longitude
  /// coordinate basis
  void fabTransformEquiangularToLonLat(
                                       const FArrayBox& a_xiFab,
                                       FArrayBox& a_rllXiFab
                                       ) const;

  /// transform a FAB of points from mapped-coordinate basis to Cartesian
  /// coordinate basis
  void fabTransformEquiangularToCartesian(
                                          const FArrayBox& a_xiFab,
                                          FArrayBox& a_xyzFab) const;

  /// transform a vector from mapped-coordinate basis to real-coordinate basis,
  /// given a point a_xi in mapped space
  void vectorTransformEquiangularToCartesian(const RealVect& a_xi,
                                             const Real * a_vecCS,
                                             Real * a_vecXYZ) const;

  /// transform a vector from real-coordinate basis to mapped-coordinate basis,
  /// given a point a_xi in mapped space
  void vectorTransformCartesianToEquiangular(const RealVect& a_xi,
                                             const Real * a_vecXYZ,
                                             Real * a_vecCS) const;

  /// transform a vector from longitude-latitude coordinate basis to
  /// mapped-coordinate basis
  void vectorTransformLatLonToEquiangular(const RealVect& a_xi,
                                          const Real * a_vecRLL,
                                          Real * a_vecCS) const;

  /// transform a FAB of vectors from latitude-longitude coordinate basis to
  /// mapped-coordinate basis
  void fabVectorTransformLatLonToEquiangular(
                                             const FArrayBox& a_xiFab,
                                             const FArrayBox& a_vecRLLFab,
                                             FArrayBox& a_vecCSFab) const;

  /// transform a FAB of vectors from mapped-coordinate basis to
  /// latitude-longitude coordinate basis
  void fabVectorTransformEquiangularToLatLon(
                                             const FArrayBox& a_xiFab,
                                             const FArrayBox& a_vecCSFab,
                                             FArrayBox& a_vecRLLFab) const;

  /// transform a FAB of SpaceDim-vectors from mapped-coordinate basis to real-coordinate basis at cell centers
  virtual void vectorTransformMappedToRealCenterFab(FArrayBox& a_vectorFab) const;

  /// transform a FAB of SpaceDim-vectors from real-coordinate basis to mapped-coordinate basis at cell centers
  virtual void vectorTransformRealToMappedCenterFab(FArrayBox& a_vectorFab) const;

  /// return row a_idir of contravariant metric
  virtual void contravariantMetric(FArrayBox& a_metric,
                                   int a_dir) const;

  /// convert components in a_csComps of vector a_csFab in coordinate-system basis to components in a_orthoComps of vector a_orthoFab in orthonormal basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void orthonormalize(
                              const FArrayBox& a_csFab,
                              FArrayBox& a_orthoFab,
                              const Box& a_box,
                              int a_idir,
                              const IntVect& a_csComps,
                              const IntVect& a_orthoComps) const;

  /// convert components in a_orthoComps of vector a_orthoFab in orthonormal basis to components in a_csComps of vector a_csFab in coordinate-system basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void deorthonormalize(
                                const FArrayBox& a_orthoFab,
                                FArrayBox& a_csFab,
                                const Box& a_box,
                                int a_idir,
                                const IntVect& a_orthoComps,
                                const IntVect& a_csComps) const;

  /// return coefficients to convert components of a vector from coordinate-system basis to orthonormal basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void getOrthonormalizingMatrix(
                                         FArrayBox& a_matrixFab,
                                         const Box& a_box,
                                         int a_idir) const;

  /// return coefficients to convert components of a vector from orthonormal basis to coordinate-system basis, at points with indices in a_box, on lines where direction a_idir is constant
  virtual void getDeorthonormalizingMatrix(
                                           FArrayBox& a_matrixFab,
                                           const Box& a_box,
                                           int a_idir) const;

  /// return radial component of curl of a_vecFab, which has SpaceDim components and is cell-centered.  Does not seem to work.
  virtual void curl(
                    const FArrayBox& a_vecFab,
                    FArrayBox& a_curlFab,
                    const Box& a_box) const;

  
  /// return radial component of curl of a_vecFab, which has SpaceDim components and is cell-centered, by converting to spherical coordinates.
  virtual void curlSpherical(
                             const FArrayBox& a_vecFab,
                             FArrayBox& a_curlFab,
                             const Box& a_box) const;

  /// access function to panel id
  int panel() const
  {
    return m_nPanel;
  }

protected:
  // Panel id (this object has different behaviour depending on the panel)
  int m_nPanel;

  // number of dimensions in real space
  int m_realDim;

  // whether to map to a flat plane; default false
  bool m_flatMap;

  // Index of the bottom-left edge of this panel in index space
  IntVect m_ix;

};

#include "NamespaceFooter.H"

#endif
