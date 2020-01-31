#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CIRCLEADVECTMULTIMAPPEDIBC_H_
#define _CIRCLEADVECTMULTIMAPPEDIBC_H_

#include  <iostream>

#include "RealVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Initial conditions for gaussian.
/**
   Parameters used in the initial conditions
     m_rr0            - initial radius
 */
class CircleAdvectMultiMappedIBC:public PhysAdvectMappedIBC
{
public:
  // type of velocity field
  enum velType
  {
    UNIFORM = 0,
    SOLIDBODY,
    NUM_VEL_TYPE
  };

  /// Null constructor
  /**
   */
  CircleAdvectMultiMappedIBC();

  /// Destructor
  /**
   */
  virtual ~CircleAdvectMultiMappedIBC()
  {
  }

  /// Sets parameters
  /**
   */
  void setParams(Real a_r0,
                 RealVect a_center,
                 RealVect a_x0=RealVect::Zero);

  /// set index of advected quantity and interval of velocity
  void setIndices(int a_indAdvected,
                  const Interval& a_velInt);

  /// set uniform velocity field
  void setUniformVel(const RealVect& a_vel);

  /// set solid-body rotation velocity field
  void setSolidBodyRotation(const RealVect& a_rotationCenter,
                            const Real a_omega);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysMappedIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_params_are_set
     set to value of m_params_are_set in the current (factory) object.
   */
  PhysMappedIBC *new_physIBC();

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_phi);

  void initializeWithJ(LevelData<FArrayBox>& a_phi);

  /// set up face-centered advection velocity
  void advVel(LevelData<FArrayBox>& a_advVel,
              Real a_time = 0.);

  /// fill ghost cell values at domain boundaries
  virtual void ghostCellBC(LevelData<FArrayBox>& a_phi);


  /// compute exact solution
  /**
     for Circle, this is a call to initialize
   */
  void exactSoln(LevelData<FArrayBox>& a_phi);

  /// Set boundary primitive values.
  /**
   */
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

protected:
  RealVect m_x0;
  RealVect m_center;
  Real m_rr0;
  bool m_params_are_set;
  int m_velType;
  RealVect m_uniformVel;
  RealVect m_rotationCenter;
  Real m_omega;

  void initializeUnified(LevelData<FArrayBox>& a_phi,
                         bool a_includeJ);

  /// Set up initial conditions
  /**
   */
  virtual void pointVal(FArrayBox& a_phi,
                        const Box& a_box,
                        const Box& a_bxBase,
                        bool a_includeJ);
};

#include "NamespaceFooter.H"

#endif
