#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TWISTEDCOORDSYS_H_
#define _TWISTEDCOORDSYS_H_

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

/// Calhoun-Leveque "twisted" Coordinate mapping (constant Jacobian)
/**
*/
class TwistedCoordSys : public FourthOrderCoordSys
{
public:
  /// default constructor
  TwistedCoordSys();

  /**
     Destructor.
   */
  virtual ~TwistedCoordSys();

  void define(const DisjointBoxLayout& a_grids,
              const ProblemDomain& a_domain,
              const RealVect& a_cellSpacing,
              const IntVect& a_ghostVect);

  void radius(const Real& a_R)
  {
    m_R = a_R;
    m_scale = m_Pi/m_R;
  }

  void twist(const Real& a_theta)
  {
    // in degrees
    m_theta = a_theta;
  }

  void regrid(const DisjointBoxLayout& a_newGrids);

  RealVect realCoord(const RealVect& a_X) const;

  RealVect mappedCoord(const RealVect& a_x) const;

  Real pointwiseJ(const RealVect& a_X) const;

  Real radius() const
  {
    return m_R;
  }

  Real twist() const
  {
    return m_theta;
  }

  void mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                            LevelData<FluxBox>& a_F);

protected:


  virtual Real dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const;

  Vector<LevelData<FluxBox>* > m_dNdXi;

  // radius for n-sphere of twisted grid
  Real m_R;

  // Amount of twist in the grid
  Real m_theta;

  // convenience
  Real m_Pi;
  Real m_scale;

  bool m_isDefined;


};

/// factory for TwistedCoordSys
/**
 */
class
TwistedCoordSysFactory : public CoordSysFactory<FArrayBox, FluxBox>
{

public:
  /// constructor
  /** baseDomain is problemDomain at coarsest level (level 0)
      vectRefRefine are refinement ratios
      baseDx is cell spacing at level 0
      maxLevel is defined by vectRefRefine.size() + 1
  */
  TwistedCoordSysFactory(const ProblemDomain& a_baseDomain,
                         const Vector<int>& a_vectRefRatios,
                         const RealVect& a_baseDx,
                         const Real& a_R,
                         const Real& a_twist);




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
  Real m_R;
  Real m_theta;
};



#endif
