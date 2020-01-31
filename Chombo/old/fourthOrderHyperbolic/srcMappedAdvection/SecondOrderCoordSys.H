#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SECONDORDERCOORDSYS_H_
#define _SECONDORDERCOORDSYS_H_

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
   The SecondOrderCoordSys class extends the CoordSys base class by
   providing API elements and implementation details specific to
   second-order discretizaitions.
   Implementations of specific coordinate systems will be provided by
   derived classes.
*/
class SecondOrderCoordSys : public CoordSys<FArrayBox, FluxBox>
{
public:
  /// default constructor
  SecondOrderCoordSys();

  /**
     Destructor.
   */
  virtual ~SecondOrderCoordSys();


  /// this define should probably be called from the derived-class define
  virtual void define(const DisjointBoxLayout& a_grids,
                      const ProblemDomain& a_domain,
                      const RealVect& a_cellSpacing,
                      const IntVect& a_ghostVect);

  virtual void regrid(const DisjointBoxLayout& a_newGrids) = 0;

  virtual RealVect realCoord(const RealVect& a_Xi) const = 0;

  virtual RealVect mappedCoord(const RealVect& a_x) const = 0;


  virtual void mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                    const LevelData<FluxBox>& a_F);


  virtual void setQuadrature(const EdgeQuadrature* a_quadrature);

  virtual const LevelData<FluxBox>& getFaceMetricTerms() const;

  virtual const LevelData<FArrayBox>& getCellVolumes() const;

  /// returns second-order averaged face-centered 1/J
  virtual const LevelData<FluxBox>& getJInverse() const;

  /// returns second-order averaged face-centered N/J
  virtual const LevelData<FluxBox>& getNJinverse() const;

  /// note that a_X is in real (not mapped) space.
  virtual Real  dXdXi(const RealVect& a_X, int a_dirX,
                      int a_dirXi) const = 0;

  /// which component directions to use when computing cell volumes
  Interval m_volInterval;

protected:

  RealVect m_dx;

  ProblemDomain m_domain;

  IntVect m_ghostVect;

  LevelData<FluxBox> m_faceMetricTerms;

  LevelData<FArrayBox> m_cellVolumes;

  LevelData<FluxBox> m_JInverse;

  LevelData<FluxBox> m_NinverseJ;

  EdgeQuadrature* m_quadraturePtr;

  /// note that a_X is in real (not mapped) space.
  virtual Real getN(const RealVect& a_X, int a_s, int a_d, int a_d1) const;

  /// note that a_X is in real (not mapped) space.
  virtual Real getNMatrixEntry(const RealVect& a_X,
                               int a_s, int a_d, int a_d1,
                               int a_row, int a_column) const;


  virtual void computeFaceFlux(LevelData<FluxBox>& a_Flux,
                               const LevelData<FluxBox>& a_F);

  virtual void defineMetricTerms();

  virtual void defineFaceMetricTerms(LevelData<FluxBox>& a_faceMetricTerms);

  virtual void incrementFaceMetricWithEdgeTerm(FArrayBox& a_faceMetrics,
                                               int a_faceDir,
                                               int a_edgeDir,
                                               const Side::LoHiSide& a_side);


  virtual void computeCellVolumes(LevelData<FArrayBox>& a_cellVolumes);

  virtual void computeJinverse(LevelData<FluxBox>& a_Jinverse);

  virtual void computeNJinverse(LevelData<FluxBox>& a_NJinverse);
};

#include "NamespaceFooter.H"
#endif
