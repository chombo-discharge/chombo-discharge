#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RTHETAZCOORDSYS_H_
#define _RTHETAZCOORDSYS_H_

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

/// Basic R-Z-Theta Coordinates
/**
   The R-Z-ThetaCoordSys class implements the FourthOrderCoordSys class
   for simple R-Z-Theta coordinates
*/
class RThetaZCoordSys : public FourthOrderCoordSys
{
public:
  /// default constructor
  RThetaZCoordSys();

  /**
     Destructor.
   */
  virtual ~RThetaZCoordSys();

  void define(const DisjointBoxLayout& a_grids,
              const ProblemDomain& a_domain,
              const RealVect& a_cellSpacing,
              const IntVect& a_ghostVect);

  void stretch(const RealVect& a_stretch)
  {
    m_stretch = a_stretch;
  }

  void regrid(const DisjointBoxLayout& a_newGrids);

  RealVect realCoord(const RealVect& a_X) const;

  RealVect mappedCoord(const RealVect& a_x) const;

  Real pointwiseJ(const RealVect& a_X) const;

  RealVect stretch() const
  {
    return m_stretch;
  }

  void mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                            LevelData<FluxBox>& a_F);

protected:


  virtual Real dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const;

  Vector<LevelData<FluxBox>* > m_dNdXi;

  // allows for stretched RThetaZ grids
  RealVect m_stretch;

  // convenience
  Real m_Pi;

  bool m_isDefined;


};

/// factory for RThetaZCoordSys
/**
 */
class
RThetaZCoordSysFactory : public CoordSysFactory<FArrayBox, FluxBox>
{

public:
  /// constructor
  /** baseDomain is problemDomain at coarsest level (level 0)
      vectRefRefine are refinement ratios
      baseDx is cell spacing at level 0
      maxLevel is defined by vectRefRefine.size() + 1
  */
  RThetaZCoordSysFactory(const ProblemDomain& a_baseDomain,
                           const Vector<int>& a_vectRefRatios,
                           const RealVect& a_baseDx,
                           const RealVect& a_stretch = RealVect::Unit);




  ///
  virtual CoordSys<FArrayBox, FluxBox>* getCoordSys(const DisjointBoxLayout& grids,
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
  RealVect m_stretch;
};


#endif
