#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERESHELLCS_H_
#define _CUBEDSPHERESHELLCS_H_

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
#include "Spline1DMapping.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

/// Multiblock coordinate system for cubed sphere shell
/// with variable-thickness spherical layers 
/**
   The CubedSphereShellCS coordinate system class is a multiblock example
   that takes place on a cubed inflated onto the surface of a sphere, and
   creates a shell with r-direction layers of variable thickness.
   NOTE: Supports only a "scaled thin atmosphere" approximation which 
   ignores the r dependencies of the horizontal grid metrics
*/
class CubedSphereShellCS : public MultiBlockCoordSys
{
public:
  /// default constructor
  CubedSphereShellCS();

  /// constructor that makes constant thickness vertical layers
  CubedSphereShellCS(int a_nResolution, int a_nLayers = 0, 
      Real a_height = 1, Real a_radius = 1);

  /**
     Destructor.
   */
  virtual ~CubedSphereShellCS();

  virtual void define(int a_nResolution, int a_nLayers = 0, 
      Real a_height = 1, Real a_radius = 1);

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  void setFlatMap(bool a_flatMap);

  void setHeight(Real a_height);

  Real getHeight() const;

  void setRadius(Real a_radius);

  Real getRadius() const;

  void setVerticalMap(RefCountedPtr<Spline1DMapping> a_map);

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

  /// convert from Lon-Lat-R coordinates to CSS computational coordinate
  /** given a_llr, find a_xi and block a_nPanel. NOTE: no spline lookup
  */
  void cssXiFromLonLat2D(RealVect a_llr,
                         RealVect& a_xi,
                         int& a_nPanel) const;

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

protected:
  // Resolution per element
  int m_nResolution;

  // Number of vertical layers (= resolution)
  int m_nLayers;

  // Height of the domain in r direction
  Real m_height;

  // Radius of the domain in r direction - considered constant
  Real m_radius;

  // Whether to map to a flat surface; default false
  bool m_flatMap;

  // Vertical mapping, should map from [0,1] to radial direction
  RefCountedPtr<Spline1DMapping> m_verticalMap;

  // define m_boundaries, called by define()
  void defineBoundaries();
};


class CubedSphereShellCSFactory : public MultiBlockCoordSysFactory
{
public:

  CubedSphereShellCSFactory()
  {
    m_flatMap = false;
  }

  virtual ~CubedSphereShellCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

  void setVerticalMap(RefCountedPtr<Spline1DMapping> a_map)
  {
    m_verticalMap = a_map;
  }

  void setHeight(Real a_height)
  {
    m_height = a_height;
  }

  void setRadius(Real a_radius)
  {
    m_radius = a_radius;
  }

  void setFlatMap(bool a_flatMap)
  {
    m_flatMap = a_flatMap;
  }

protected:

  Real m_height;
  Real m_radius;
  bool m_flatMap;
  RefCountedPtr<Spline1DMapping> m_verticalMap;
};

#include "NamespaceFooter.H"

#endif
