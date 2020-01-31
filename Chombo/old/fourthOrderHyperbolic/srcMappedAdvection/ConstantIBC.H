#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTANTIBC_H_
#define _CONSTANTIBC_H_

#include  <iostream>

#include "RealVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BasicIBC.H"
#include "UsingNamespace.H"

/// Initial conditions for constant.
/**

 */
class ConstantIBC:public BasicIBC
{
public:
  // type of velocity field
  enum  velType
  {
    UNIFORM = 0,
    TRANSOSC,
    SOLIDBODY,
    NUM_VEL_TYPE
  };

  /// Null constructor
  /**
   */
  ConstantIBC();

  /// Destructor
  /**
   */
  ~ConstantIBC()
  {
  }

  /// Sets parameters
  /**
   */
  void setParams(Real a_mag);

  /// set uniform velocity field
  void setUniformVel(const RealVect& a_vel);

  /// set solid-body rotation velocity field
  void setSolidBodyRotation(const RealVect& a_rotationCenter,
                            const Real a_omega);

  /// set parameters for translating-oscilating velocity
  void setTranslatingOscillation(const RealVect& a_translation_vel,
                                 const Real a_oscillation_amplitude);

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
              Real a_dXi, Real a_time);

  /// fill ghost cell values at domain boundaries
  virtual void ghostCellBC(LevelData<FArrayBox>& a_phi,
                           const ProblemDomain& a_domain,
                           const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                           Real a_dx,
                           Real a_time);


  /// compute exact solution
  /**
   */
  void exactSoln(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain,
                 const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                 Real a_dx,
                 Real a_time=0.0);


protected:
  Real m_mag;

  int m_velType;
  RealVect m_uniformVel;
  RealVect m_rotationCenter;
  Real m_omega;
  Real m_oscAmp;

  LevelData<FArrayBox> m_fineJ;
  LevelData<FArrayBox> m_crseJ;

  bool m_params_are_set;

  /// Set up initial conditions
  /**
   */
  void pointVal(FArrayBox& a_phi,
                const ProblemDomain& a_domain,
                const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                const Box& a_box,
                bool a_includeJ,
                Real a_dx,
                Real time);

};

#endif
