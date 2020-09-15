#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// GaussianAdvectFun.H
// petermc, 4 Feb 2011

#ifndef _GAUSSIANADVECTFUN_H_
#define _GAUSSIANADVECTFUN_H_

#include "FunctionOnSpace.H"
// #include "GaussianAdvectMultiMappedIBC.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class GaussianAdvectFun : public FunctionOnSpace
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves GaussianAdvectFun undefined.
  */
  GaussianAdvectFun();

  ///
  /**
     Full constructor.  Assumes solid body.
  */
  GaussianAdvectFun(Real a_r0,
                    RealVect a_center,
                    RealVect a_x0,
                    Real a_initialTime,
                    RealVect a_rotationCenter,
                    Real a_omega)
  {
    define(a_r0, a_center, a_x0, a_initialTime, a_rotationCenter, a_omega);
  }

  ///
  /** Destructor.
   */
  virtual ~GaussianAdvectFun();

  virtual void define(Real a_r0,
                      RealVect a_center,
                      RealVect a_x0,
                      Real a_initialTime,
                      RealVect a_rotationCenter,
                      Real a_omega);


  void setCoordSys(MultiBlockCoordSys* a_coordSysPtr);

  ///
  /** Set a_funFab on a_bx to function values at points with
      physical coordinates a_coordsFab.
   */
  virtual void setFunctionFromPhysical(/// function values, 1 component
                                       FArrayBox&         a_funFab,
                                       /// indices on which to compute function
                                       const Box&         a_bx,
                                       /// physical coordinates, SpaceDim components
                                       const FArrayBox&   a_coordsFab);

  ///
  /** Return maximum value of function on the multiblock domain.
   */
  virtual Real setFunctionMax(/// physical side length of central cube
                              Real   a_bxWidth,
                              /// outer radius of circle
                              Real   a_outerRadius);

protected:

  //  GaussianAdvectMultiMappedIBC* m_ibc;
  Real m_r0;
  Real m_rr0powInv2;
  Real m_rr0powInvD;
  RealVect m_center;
  RealVect m_x0;
  Real m_initialTime;
  RealVect m_rotationCenter;
  Real m_omega;
  MultiBlockCoordSys* m_coordSysPtr;

private:

  void operator = (const GaussianAdvectFun&)
  {
  }

  GaussianAdvectFun(const GaussianAdvectFun&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
