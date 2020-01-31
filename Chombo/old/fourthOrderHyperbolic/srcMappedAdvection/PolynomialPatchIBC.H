#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POLYNOMIALPATCHIBC_H_
#define _POLYNOMIALPATCHIBC_H_

#include  <iostream>

#include "RealVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BasicIBC.H"
#include "UsingNamespace.H"


/// Initial conditions for polynomialPatch.
/**
   Parameters used in the initial conditions
   m_mag,
   m_radius
   m_center
 */
class PolynomialPatchIBC:public BasicIBC
{
public:
  // type of velocity field
  enum  velType
  {
    uniform = 0,
    NUM_VEL_TYPE
  };

  /// Null constructor
  /**
   */
  PolynomialPatchIBC();

  /// Destructor
  /**
   */
  ~PolynomialPatchIBC()
  {
  }

  /// Sets parameters
  /**
   */
  void setParams(Real a_rad, RealVect a_center,
                 Real a_mag);

  /// set uniform velocity field
  void setUniformVel(const RealVect& a_vel);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_params_are_set
     set to value of m_params_are_set in the current (factory) object.
   */
  BasicIBC *new_basicIBC();

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_phi,
                  const ProblemDomain& a_domain,
                  const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                  Real a_dx,
                  Real time=0.0);

  /// set up face-centered advection velocity
  void advVel(LevelData<FArrayBox>& a_advVel,
              const ProblemDomain& a_domain,
              const CoordSys<FArrayBox,FluxBox>& a_coordSys,
              Real a_dx,
              Real a_time);

  /// fill ghost cell values at domain boundaries
  virtual void ghostCellBC(LevelData<FArrayBox>& a_phi,
                           const ProblemDomain& a_domain,
                           const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                           Real a_dx,
                           Real a_time);


  /// compute exact solution
  /**
     for PolynomialPatch, we can compute exact solution analytically
   */
  void exactSoln(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain,
                 const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                 Real a_dx,
                 Real a_time=0.0);


protected:
  Real m_rad;
  RealVect m_center;
  Real m_mag;
  bool m_params_are_set;
  int m_velType;
  RealVect m_uniformVel;

  // polynomial coefficients
  Real m_A, m_B, m_C, m_D, m_E, m_F;

  /// solution at a given time
  /**
   */
  virtual void pointVal(FArrayBox& a_phi,
                        const ProblemDomain& a_domain,
                        const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                        const Box& a_box,
                        bool a_includeJ,
                        Real a_dx,
                        Real a_time);

};

#endif
