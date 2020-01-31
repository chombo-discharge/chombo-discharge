#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERCOORDSYS_H_
#define _FOURTHORDERCOORDSYS_H_

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "RealVect.H"
#include "IntVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "CoordSys.H"
#include "EdgeQuadrature.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

/// Virtual base class encapsulating fourth-order CoordSys
/**
   The FourthOrderCoordSys class extends the CoordSys base class by
   providing API elements and implementation details specific to
   fourth-order discretizaitions.
   Implementations of specific coordinate systems will be provided by
   derived classes.
*/
class FourthOrderCoordSys : public CoordSys<FArrayBox, FluxBox>
{
public:
  /// default constructor
  FourthOrderCoordSys();

  /**
     Destructor.
   */
  virtual ~FourthOrderCoordSys();


  /// this define should probably be called from the derived-class define
  virtual void define(const DisjointBoxLayout& a_grids,
                      const ProblemDomain& a_domain,
                      const RealVect& a_cellSpacing,
                      const IntVect& a_ghostVect);

  virtual void regrid(const DisjointBoxLayout& a_newGrids) = 0;

  virtual RealVect realCoord(const RealVect& a_Xi) const = 0;

  virtual RealVect mappedCoord(const RealVect& a_x) const = 0;

  // default implementation just uses
  virtual void realCoord(FArrayBox& a_x, const FArrayBox& a_Xi) const;

  virtual void mappedCoord(FArrayBox& a_Xi, const FArrayBox& a_x) const;


  virtual Real pointwiseJ(const RealVect& a_X) const;

  /// Jacobian evaluated at index locations of a_J
  virtual void pointwiseJ(FArrayBox& a_J) const;

  virtual void pointwiseJ(FArrayBox& a_J, const FArrayBox& a_Xi) const;

  virtual void mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                    const LevelData<FluxBox>& a_F) const;

  virtual void simpleDivergence(LevelData<FArrayBox>& a_divF,
                                const LevelData<FluxBox>& a_F) const;

  virtual void setQuadrature(const EdgeQuadrature* a_quadrature);

  /// reuturns reference to face-averaged N
  virtual const LevelData<FluxBox>& getFaceMetricTerms() const;

  /// returns reference to tangential gradients of N
  virtual const LevelData<FluxBox>& getTanGranNTerms() const
  {
    return m_tanGradNTerms;
  }

  virtual int tanGradComp(const int a_faceDir, const int a_tanDir,
                          const int a_comp) const;

  virtual const LevelData<FArrayBox>& getCellVolumes() const;

  /// returns fourth-order averaged cell-centered J
  virtual const LevelData<FArrayBox>& getJ() const;

  /// returns fourth-order averaged face-centered 1/J
  virtual const LevelData<FluxBox>& getJInverse() const;

  /// returns fourth-order averaged face-centered N/J
  virtual const LevelData<FluxBox>& getNJinverse() const;


  virtual void computeTangentialGrad(LevelData<FluxBox>& a_gradPhi,
                                     const LevelData<FluxBox>& a_phiFace) const;


  /// helper function
  /**  given face-centered phi, increment with Laplacian(phi) correction
       in place to compute 4th-order face averages of phi
  */
//  virtual void fourthOrderAverage(LevelData<FluxBox>& a_phi) const;


  /// computes 4th-order average of product = F*G
  virtual void computeProductAverage(LevelData<FluxBox>& a_product,
                                     const LevelData<FluxBox>& a_F,
                                     const LevelData<FluxBox>& a_G,
                                     bool a_fourthOrder = true) const;


  /// computes 4th-order average of product = N^T*F
  /** if a_fourthOrder is false, then only do a second-order dot product
   */
  virtual void computeMetricTermProductAverage(LevelData<FluxBox>& a_product,
                                               const LevelData<FluxBox>& a_F,
                                               bool a_fourthOrder=true) const;

  /// index function into face-centered metric terms
  /** returns which component of the face-centered metric term in
      which N^row_col is stored
  */
  static inline int getMetricTermComponent(const int a_row, const int a_col)
  {
    return (SpaceDim*a_col) + a_row;
  }

  /// note that a_X is in real (not mapped) space.
  virtual Real  dXdXi(const RealVect& a_X, int a_dirX,
                      int a_dirXi) const = 0;

  /// note that a_X is in real (not mapped) space.
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
                     const FArrayBox& a_X,
                     int a_destComp,
                     int a_dirX,
                     int a_dirXi) const;

  /// which component directions to use when computing cell volumes
  Interval m_volInterval;

  /// access function to simplify things -- returns mapped-space cell spacing
  const RealVect& dx() const
  {
    return m_dx;
  }

  /// access function, returns mapped-space problem domain
  const ProblemDomain& domain() const
  {
    return m_domain;
  }

  /// access function -- returns const reference to grids
  const DisjointBoxLayout& grids() const
  {
    return m_grids;
  }

  /// access function -- returns const reference to number of ghost cells
  const IntVect& ghostVect() const
  {
    return m_ghostVect;
  }

protected:

  DisjointBoxLayout m_grids;

  RealVect m_dx;

  ProblemDomain m_domain;

  IntVect m_ghostVect;

  LevelData<FluxBox> m_faceMetricTerms;

  LevelData<FluxBox> m_tanGradNTerms;

  LevelData<FArrayBox> m_cellVolumes;

  LevelData<FArrayBox> m_J;

  LevelData<FluxBox> m_JInverse;

  LevelData<FluxBox> m_NinverseJ;

  EdgeQuadrature* m_quadraturePtr;

  /// note that a_X is in real (not mapped) space.
  virtual Real getN(const RealVect& a_X, int a_s, int a_d, int a_d1) const;

  /// note that a_X is in real (not mapped) space.
  virtual Real getNMatrixEntry(const RealVect& a_X,
                               int a_s, int a_d, int a_d1,
                               int a_row, int a_column) const;


  virtual void defineFaceMetricTerms(LevelData<FluxBox>& a_faceMetricTerms);

  virtual void incrementFaceMetricWithEdgeTerm(FArrayBox& a_faceMetrics,
                                               int a_faceDir,
                                               int a_edgeDir,
                                               const Side::LoHiSide& a_side);

  virtual void computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                           int a_faceDir,
                                           int a_dDir);

  virtual void computeCellVolumes(LevelData<FArrayBox>& a_cellVolumes);

  virtual void computeJ(LevelData<FArrayBox>& a_J);

  virtual void computeJinverse(LevelData<FluxBox>& a_Jinverse);

  virtual void computeNJinverse(LevelData<FluxBox>& a_NJinverse);
};

#include "NamespaceFooter.H"

#endif

