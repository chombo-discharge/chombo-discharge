#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WARPEDCOORDSYS_H_
#define _WARPEDCOORDSYS_H_

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
#include "FourthOrderCoordSys.H"

#include "UsingNamespace.H"

/// "warped" Coordinate mapping (non-constant Jacobian)
/**
*/
class WarpedCoordSys : public FourthOrderCoordSys
{
public:
  /// default constructor
  WarpedCoordSys();

  /**
     Destructor.
   */
  virtual ~WarpedCoordSys();

  void define(const DisjointBoxLayout& a_grids,
              const ProblemDomain& a_domain,
              const RealVect& a_cellSpacing,
              const IntVect& a_ghostVect);

  void scale(const RealVect& a_scale);
  void relative_tolerance(const Real a_rtol);
  void absolute_tolerance(const Real a_atol);
  void maximum_iterations(const int a_imax);

  void regrid(const DisjointBoxLayout& a_newGrids);

  RealVect realCoord(const RealVect& a_X) const;

  RealVect mappedCoord(const RealVect& a_x) const;

  Real pointwiseJ(const RealVect& a_X) const;

  RealVect scale() const
  {
    return m_scale;
  }

  Real relative_tolerance() const
  {
    return m_rtol;
  }

  Real absolute_tolerance() const
  {
    return m_atol;
  }

  int maximum_iterations() const
  {
    return m_imax;
  }

  void mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                            LevelData<FluxBox>& a_F);

protected:

  // RHS of fixed-point iteration for inverse mapping
  Real g( const Real xi, const RealVect& x ) const;

  virtual Real dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const;

  Vector<LevelData<FluxBox>* > m_dNdXi;

  // scale for displacement
  RealVect m_scale;

  // parameters for fixed-point iteration for inverse mapping
  Real m_rtol;
  Real m_atol;
  int m_imax;

  // convenience
  Real m_twoPi;
  int m_rootDir;

  bool m_isDefined;
};

/// factory for WarpedCoordSys
/**
 */
#ifdef CH_USE_FLOAT
#  define RTOL 1.e-6
#  define ATOL 1.e-7
#else
#  define RTOL 1.e-12
#  define ATOL 1.e-13
#endif


class
WarpedCoordSysFactory : public CoordSysFactory<FArrayBox, FluxBox>
{

public:
  /// constructor
  /** baseDomain is problemDomain at coarsest level (level 0)
      vectRefRefine are refinement ratios
      baseDx is cell spacing at level 0
      maxLevel is defined by vectRefRefine.size() + 1
  */
  WarpedCoordSysFactory(const ProblemDomain& a_baseDomain,
                        const Vector<int>& a_vectRefRatios,
                        const RealVect& a_baseDx,
                        const RealVect& a_scale,
                        const Real a_rtol=RTOL,
                        const Real a_atol=ATOL,
                        const int a_imax=100);




  ///
  virtual CoordSys<FArrayBox, FluxBox>* getCoordSys(
     const DisjointBoxLayout& grids,
     const ProblemDomain& levelDomain,
     const IntVect& a_ghostVect) const;

  /// interval of which components over which to compute cell volumes
  Interval m_volInterval;

protected:

  ///
  Vector<ProblemDomain> m_vectDomains;

  ///
  Vector<int> m_vectRefRatios;

  ///
  Vector<RealVect> m_dxVect;

  ///
  RealVect m_scale;
  Real m_rtol;
  Real m_atol;
  int m_imax;
};



#endif
