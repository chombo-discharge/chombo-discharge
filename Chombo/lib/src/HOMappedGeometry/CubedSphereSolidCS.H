#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERESOLIDCS_H_
#define _CUBEDSPHERESOLIDCS_H_

//#include <utility>
// using std::pair;

#include "REAL.H"
#include "RealVect.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Multiblock coordinate system for cubed sphere shell
/// with variable-thickness spherical layers 
/**
   The CubedSphereSolidCS coordinate system class is a multiblock example
   that takes place on a 3D spherical region.
   There are an inner cube and six sectors, each sector centered on
   a coordinate axis.

   From: Diener, Dorband, Schnetter and Tiglio,
   "Optimized High-Order Derivative and Dissipation Operators
   Satisfying Summation by Parts, and Applications in Three-dimensional
   Multi-block Evolutions", J. Sci. Comp. 32:1 (2007)
   DOI: 10.1007/s10915-006-9123-7
*/
class CubedSphereSolidCS : public MultiBlockCoordSys
{
public:

  /// default constructor
  CubedSphereSolidCS();

  /// constructor with resolution
  CubedSphereSolidCS(int a_nResolution, Real a_r0 = 1., Real a_r1 = 3.);

  /**
     Destructor.
   */
  virtual ~CubedSphereSolidCS();

  virtual void define(int a_nResolution, Real a_r0 = 1., Real a_r1 = 3.);

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  void setR0(Real a_r0);

  void setR1(Real a_r1);

  Real getR0() const { return m_r0; }

  Real getR1() const { return m_r1; }

  virtual void blockRemapping(RealVect& a_xi_valid,
                              int& a_n_valid,
                              const RealVect& a_xiSrc,
                              int a_nSrc) const;

protected:

  const static int NUMBLOCKS = 7;

  const static int XPOS = 0;
  const static int YPOS = 1;
  const static int XNEG = 2;
  const static int YNEG = 3;
  const static int ZPOS = 4;
  const static int ZNEG = 5;
  const static int CENTRAL = 6;

  // Resolution per element
  int m_nResolution;

  // half-width of inner cube in physical space
  Real m_r0;

  // radius of sphere in physical space
  Real m_r1;

  // define m_boundaries, called by define()
  void defineBoundaries();
};


class CubedSphereSolidCSFactory : public MultiBlockCoordSysFactory
{
public:

  CubedSphereSolidCSFactory()
  {
    m_r0 = 1.;
    m_r1 = 3.;
  }

  virtual ~CubedSphereSolidCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

  void setR0(Real a_r0)
  {
    m_r0 = a_r0;
  }

  void setR1(Real a_r1)
  {
    m_r1 = a_r1;
  }

protected:

  Real m_r0;

  Real m_r1;
};

#include "NamespaceFooter.H"

#endif
