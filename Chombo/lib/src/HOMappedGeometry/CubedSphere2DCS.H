#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERE2DCS_H_
#define _CUBEDSPHERE2DCS_H_

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
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

/// simple multiblock coordinate system for a cylinder
/**
   The CubedSphere2DCS coordinate system class is a multiblock example
   that takes place on a cubed inflated onto the surface of a sphere.
*/
class CubedSphere2DCS : public MultiBlockCoordSys
{
public:
  /// default constructor
  CubedSphere2DCS();

  CubedSphere2DCS(int a_nResolution,
                  int a_nLayers = 0);

  /**
     Destructor.
   */
  virtual ~CubedSphere2DCS();

  virtual void define(int a_nResolution,
                      int a_nLayers = 0);

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  void setFlatMap(bool a_flatMap);

  void setRadii(Real a_rMin,
                Real a_rMax);

  void regrid(const DisjointBoxLayout& a_newGrids);

  RealVect realCoord(const RealVect& a_X) const;

  RealVect mappedCoord(const RealVect& a_x) const;

  /// Calculate the pointwise transformation matrix using gnomonic
  /// coordinates on the source panel (a_dX, a_dY).
  void vecTransMatrix_Pt(
    int a_nSourcePanel,
    int a_nDestPanel,
    Real a_dX,
    Real a_dY,
    Real* a_dT
  ) const;

  // Calculate the element averaged transformation matrix, and gradient
  // of the transformation matrix
  void vecTransMatrix(
        int a_nSourcePanel,
    int a_nDestPanel,
    int a_nAlphaIx,
    int a_nBetaIx,
    Real a_dDeltaA,
    Real* a_dAvgT,
    Real* a_dAlphaGradT,
    Real* a_dBetaGradT
  ) const;

  // Calculate the pointwise vector transformation
  VectorTransformation vectorBlockTransformation(
    int a_nDst,
    const RealVect& a_xiSrc,
    int a_nSrc
  ) const;

  /// convert from mapped coordinate to Cartesian coordinates
  /** given a_xiSrc in block a_nSrc, return 3D Cartesian coordinates of point.
  */
  void cartesianFromMapped(Real& xx,
                           Real& yy,
                           Real& zz,
                           const RealVect& a_xiSrc,
                           int a_nSrc) const;

  /// convert from Cartesian coordinates to mapped coordinate
  /** given 3D Cartesian coordinates of point and destination block
      returns mapped coordinates of that point.  Radius information is
      discarded.
  */
  void mappedFromCartesian(RealVect& a_xiDst,
                           int a_nDst,
                           Real xx,
                           Real yy,
                           Real zz) const;

  /// block mapping conversion function
  /** given a_xiSrc in block a_nSrc, returns which block (a_n_valid)
      contains that location in its valid region and its
      location (a_xi_valid) in that block's computational coordinate space.
  */
  virtual void blockRemapping(RealVect& a_xi_valid,
                              int& a_n_valid,
                              const RealVect& a_xiSrc,
                              int a_nSrc) const;

  // general block mapping conversion function
  /** given mapped coordinate a_xiSrc of a point in block a_nSrc return
      mapped coordinates in block a_nDst.
  */
  virtual RealVect blockRemappingGeneral(int a_nDst,
                                         const RealVect & a_xiSrc,
                                         int a_nSrc) const;

  virtual Vector<RealVect> displacements(const Vector<RealVect>&   a_dstCoords,
                                         const Vector<int>&        a_dstBlocks,
                                         const RealVect&           a_srcCoords,
                                         int                       a_srcBlock) const;

  /// Have separate volume flux on each block.  Normally does nothing, but cubed-sphere coordinate system needs it.
  virtual void separateVolFlux(LevelData<FluxBox>& a_flux) const;

  /// Function to convert from 0:1 to m_rMin:m_rMax.
  virtual Real radialMappedToReal(Real a_rMapped) const;

  /// Function to convert from m_rMin:m_rMax to 0:1.
  virtual Real radialRealToMapped(Real a_rReal) const;

protected:
  // Resolution per element
  int m_nResolution;

  // Number of vertical layers (= resolution)
  int m_nLayers;

  // Whether to map to a flat surface; default false
  bool m_flatMap;

  // These are used only in 3D.
  Real m_rMin, m_rMax;

  // Used only in 3D.  Call only after m_rMin and m_rMax have been set.
  void setFortranCommon();

  // define m_boundaries, called by define()
  void defineBoundaries();
};


class CubedSphere2DCSFactory : public MultiBlockCoordSysFactory
{
public:

  CubedSphere2DCSFactory()
  {
    m_flatMap = false;
  }

  virtual ~CubedSphere2DCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

  void setFlatMap(bool a_flatMap)
  {
    m_flatMap = a_flatMap;
  }

  void setRadii(Real a_rMin,
                Real a_rMax)
  {
    m_rMin = a_rMin;
    m_rMax = a_rMax;
  }

protected:

  bool m_flatMap;

  // These are used only in 3D.
  Real m_rMin, m_rMax;
};

#include "NamespaceFooter.H"

#endif
